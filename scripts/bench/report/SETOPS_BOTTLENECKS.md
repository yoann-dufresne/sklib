# Goulots d'étranglement temps des requêtes ensemblistes — analyse détaillée

**Date :** 2026-06-06 · **Machine :** CPU ~4.5 GHz (core), Linux, 1 cœur épinglé (`taskset -c 0`) ·
**Build :** Release `-O3 -march=native` (timing) + un build `Profile -fno-inline` (callgrind) ·
**k = 21, m = 11, buckets = 4096** (sauf sweeps) · **Données : génomes réels** (E. coli K12/Sakai, levure,
*C. elegans*, humain chr21/chr22/chr1) + copies mutées contrôlées pour piloter l'overlap.

> Ce rapport **profile l'intérieur de sklib** pour localiser le temps. Il complète `SETOPS_REPORT.md`
> (comparaison sklib vs KMC vs CBL). Toutes les mesures sont **cache chaud** (voir §I/O). Méthodes :
> décomposition macro (`_size` vs matérialisé), `perf` (cycles + instructions, attribution par fonction),
> `callgrind` (graphe d'appels exact, `-fno-inline`), et sweeps de paramètres (m, buckets, k).

---

## TL;DR — où passe le temps

1. **Pour `intersection`/`union` matérialisées, 85–91 % du temps est le POST-MERGE** : reconstruire des
   super-k-mers à partir des k-mers du résultat (`generate_sorted_list_from_enumeration`, appelé par bucket).
   Le merge lui-même (lecture + scan + comparaison) n'est que **10–15 %**.
2. **Le temps total ∝ taille de la SORTIE** ; le merge ∝ taille de l'ENTRÉE. Contrôle décisif : `diff(A,A)`
   (sortie vide) à chr1 = **2.98 s ≈ le merge seul `_size` = 2.94 s** (matérialisation sans sortie ≈ merge).
3. **Dans le post-merge, 3 fonctions concentrent ~70 % du temps total** : `colinear_chaining` ~30 %,
   `get_candidate_overlaps` ~27 %, `sort_column` ~14 % — la machinerie d'assemblage des super-k-mers,
   ré-exécutée intégralement sur la sortie, **alors que le merge l'a déjà produite triée par colonne**.
4. **Pour la cardinalité seule (`_size`)**, tout le coût est le merge : ~52 % `kmer_compare`,
   ~31 % hole-scan `set_next_valid`/`has_valid_kmer` (la boucle qui re-balaie chaque colonne), ~11 % driver.
5. **Leviers de réglage faibles** : m (merge −37 % mais matérialisation plate), buckets (−19 % à 16384),
   k (plat, +15 % seulement à k=41 → uint128). I/O ≈ 0 à chaud. Alloc < 1 %.
6. **Les buckets sont indépendants et le plus gros ne pèse que 0.2–0.3 % du travail** → un parallélisme
   par bucket **à ordonnancement dynamique** scalerait quasi-linéairement.

---

## 1. Décomposition macro : `_size` (merge) vs matérialisé (merge + post-merge)

`OP_size` ne matérialise rien (sink compteur) ; `OP` matérialise (collecte + re-compaction + écriture).
`delta = t_mat − t_size = collect(get_skmer_of_kmer) + re-compaction + write`. Médiane de 5–7 runs, cache chaud.
Données : `report/data/bottleneck_decompose.csv`.

| paire | op | sortie (k) | t_size (s) | t_mat (s) | **part post-merge** |
|---|---|--:|--:|--:|--:|
| ecoli 4.5M vs mut1% | inter | 3.69M | 0.080 | 0.568 | **0.86** |
| ecoli | union | 5.42M | 0.080 | 0.828 | **0.90** |
| ecoli | diff | 0.85M | 0.081 | 0.196 | 0.59 |
| celegans 91M vs mut1% | inter | 74.4M | 2.06 | 15.85 | **0.87** |
| celegans | union | 109.9M | 2.04 | 22.02 | **0.91** |
| chr1 190M **self** | inter | 189.8M | 2.94 | 34.40 | **0.91** |
| chr1 190M **self** | **diff (sortie vide)** | **0** | 2.94 | **2.98** | **0.01** |

**Lecture.** Le `_size` (merge) est **quasi constant selon l'opération** (ecoli : 0.080/0.080/0.081 pour
∩/∪/diff) — le merge fait le même travail quelle que soit l'op. La matérialisation, elle, suit la **sortie** :
sortie nulle ⇒ post-merge nul (diff self) ; grosse sortie (union) ⇒ post-merge maximal.

### Overlap (taille constante, ecoli) — la sortie pilote tout

| B | Jaccard ∩ | ∩ sortie | part post-merge ∩ | diff sortie | part post-merge diff |
|---|--:|--:|--:|--:|--:|
| self (100 %) | 4.54M | 4.54M | 0.91 | 0 | **0.03** |
| mut 0.5 % | 4.10M | 4.10M | 0.89 | 0.44M | 0.46 |
| mut 1 % | 3.69M | 3.69M | 0.86 | 0.85M | 0.59 |
| mut 5 % | 1.56M | 1.56M | 0.65 | 2.98M | 0.78 |
| Sakai (réel) | 3.06M | 3.06M | 0.79 | 1.48M | 0.65 |

Quand l'overlap baisse, la sortie de ∩ chute et **sa part post-merge chute avec elle** (0.91 → 0.65) ; à
l'inverse la sortie de diff grossit et sa part post-merge monte. La part post-merge **est** la part de sortie.

---

## 2. Scaling & débit (4.5M → 190M k-mers, ×42)

| génome | \|A\| Mk | merge (s) | matér. (s) | sortie Mk | **merge ns/entrée** | **post-merge ns/sortie** |
|---|--:|--:|--:|--:|--:|--:|
| ecoli | 4.5 | 0.080 | 0.568 | 3.7 | 17.5 | 132 |
| yeast | 11.5 | 0.218 | 1.462 | 9.4 | 18.9 | 133 |
| chr21 | 32.7 | 0.699 | 4.608 | 26.9 | 21.4 | 146 |
| celegans | 91.2 | 2.063 | 15.85 | 74.4 | 22.6 | 185 |
| chr1 (self) | 189.8 | 2.944 | 34.40 | 189.8 | 15.5 | 166 |

- **Merge ≈ linéaire** en entrée (~15–23 ns/k-mer ; légère hausse avec la taille = effets cache).
- **Post-merge ≈ linéaire** en sortie mais **~7–10× plus cher par élément** (~132–185 ns/k-mer), et son coût
  unitaire **croît ~40 %** d'ecoli à celegans → pression cache quand les buckets grossissent.

---

## 3. Attribution par fonction (`perf`, build Release, celegans ∩ matérialisé, ~30 k échantillons)

Agrégateur : `scripts/bench/bottleneck/perf_categorize.py`. Brut : `report/data/bottleneck_perf_summary.txt`.

| Phase / fonction | ∩ matér. | ∪ matér. | chr1 ∩ | ∩ `_size` |
|---|--:|--:|--:|--:|
| **RE-COMPACTION (rebuild super-k-mers)** | **81.1 %** | **85.8 %** | **88.4 %** | — |
| `colinear_chaining` (chaînage colinéaire, Fenwick DP) | 29.3 % | 31.2 % | 32.1 % | — |
| `get_candidate_overlaps` (détection d'overlaps) | 26.4 % | 28.0 % | 28.8 % | — |
| `sort_column` (tri par colonne) | 13.5 % | 14.1 % | 14.8 % | — |
| `merge_LList_column` (fusion des colonnes) | 5.6 % | 5.9 % | — | — |
| `Virtual_skmer` / `std::sort` / vecteurs | 6.0 % | 6.3 % | — | — |
| **MERGE (scan + compare + collect)** | 14.7 % | 10.5 % | 8.6 % | **94.6 %** |
| I/O (cache chaud) | ~0 % | ~0 % | ~0 % | ~0 % |
| alloc / mem | 0.5 % | 0.3 % | — | — |

`colinear_chaining` + `get_candidate_overlaps` + `sort_column` = **~70 % du temps total**, à toutes les échelles.

### 3-bis. Cycles vs instructions (même binaire Release) — qui *stalle* ?

Une fonction dont la part de **cycles** (temps) dépasse sa part d'**instructions** (travail) est limitée par
la mémoire/latence. IPC global = 2.06 (pas globalement memory-bound).

| fonction | cycles % | instr % | **cyc/inst** | interprétation |
|---|--:|--:|--:|---|
| `colinear_chaining` | 29.3 | 24.8 | **1.18×** | travail réel **+** stall modéré (Fenwick + tri + compression de coord. → accès dispersés) |
| `get_candidate_overlaps` | 26.4 | 23.9 | 1.10× | surtout du travail (extract + tri + `lower_bound`) |
| merge driver (+collect/`get_skmer_of_kmer`) | 14.7 | 9.3 | **1.58×** | le plus stallé : sauts de curseur + écriture dispersée des skmers collectés |
| `sort_column` | 13.5 | 25.3 | **0.53×** | **dense en instructions mais file vite** (std::sort bien optimisé) → pas un goulot proportionnel à son volume |
| `merge_LList_column` | 5.6 | 6.5 | 0.86× | — |

> Leçon : ne pas se fier au comptage d'instructions (callgrind) pour le *temps*. Ex. : à `-fno-inline`,
> `colinear_chaining` semblait ne peser que 6 % des instructions — artefact de build ; sur le vrai binaire
> Release il pèse 29 % du temps. Et `sort_column`, lourd en instructions, est en réalité rapide (0.53 cyc/inst).

---

## 4. Anatomie du merge (`callgrind`, `_size`, niveau instruction)

Pour la cardinalité seule, **tout** le coût est le merge. Dissection (instructions, `-fno-inline` ; à pondérer
par §3-bis — les pas de scan sont cache-résidents, donc *moins* chers en temps qu'en instructions) :

| composant | part instr. | rôle |
|---|--:|---|
| `set_next_valid` + `has_valid_kmer` (**hole-scan**) | **~31 %** | `merge_columns` **re-balaie tout le span depuis 0 à chaque colonne** (k−m+1 = 11 passes), en sautant les records sans k-mer à cette colonne |
| `kmer_compare` (+ `pair` copie, `&=` ×2, `<`, `>`) | **~52 %** | comparaison des k-mers sur le front de fusion ; copie le `pair` et masque 2× par appel |
| driver `merge_columns` + sink | ~11 % | boucle de fusion à deux curseurs |

Le m-sweep (§5) confirme **en temps réel** que le hole-scan coûte : le merge baisse de **37 %** quand on passe
de 12 colonnes (m=9) à 2 colonnes (m=19).

---

## 5. Sweeps de paramètres — ce qui n'aide (presque) pas

### m (k=21, levure ; `bottleneck_sweep_m.csv`)
| m | colonnes | records | t_size | t_mat | part post-merge |
|--:|--:|--:|--:|--:|--:|
| 9 | 12 | 2.91M | 0.260 | 1.593 | 0.84 |
| 11 | 10 | 2.46M | 0.215 | 1.450 | 0.85 |
| 13 | 8 | 2.43M | 0.179 | 1.438 | 0.88 |
| 19 | 2 | 5.75M | 0.163 | 1.342 | 0.88 |

Le **merge** baisse avec m (moins de colonnes à scanner) **mais la matérialisation reste plate** (moins de
colonnes ⇄ plus de records, courts). **m n'est pas un levier** pour les set-ops.

### buckets (k=21 m=11, levure ; `bottleneck_sweep_buckets.csv`)
| buckets | t_size | t_mat | part post-merge |
|--:|--:|--:|--:|
| 64 | 0.258 | 1.898 | 0.86 |
| 1024 | 0.255 | 1.676 | 0.85 |
| 4096 (défaut) | 0.235 | 1.601 | 0.85 |
| 16384 | 0.229 | **1.541** | 0.85 |
| 65536 | 0.249 | 1.544 | 0.84 |

**Levier modéré** : +buckets ⇒ tris/chaînage par bucket plus petits (effet n·log(n/B)), **−19 %** de 64 à
16384, puis rendements décroissants (surcoût fixe par bucket + bits de quotient). **16384 ≈ optimum** ici.

### k (ecoli, m=11 ; `bottleneck_sweep_k.csv`) — effet largeur de record
| k | largeur store | octets/rec | records | t_size | t_mat |
|--:|--:|--:|--:|--:|--:|
| 15 | uint32 (4) | 16 | 1.72M | 0.072 | 0.531 |
| 21 | uint32 (4) | 16 | 0.88M | 0.074 | 0.519 |
| 31 | uint64 (8) | 24 | 0.48M | 0.072 | 0.508 |
| 41 | **uint128 (16)** | **48** | 0.33M | 0.077 | **0.593** |

Quasi plat (k↑ ⇒ records moins nombreux mais plus larges, ça se compense), **+15 % seulement à k=41** quand
le record bascule en uint128 (256 bits, comparaisons plus chères).

---

## 6. Équilibre des buckets ⇒ faisabilité du parallélisme

| liste | buckets | non vides | médiane | moyenne | max | **max/moy** | **max/total** |
|---|--:|--:|--:|--:|--:|--:|--:|
| ecoli | 4096 | 2525 | 6 | 215 | 2037 | 9.5× | **0.23 %** |
| celegans | 4096 | 2924 | 129 | 6389 | 87 891 | 13.8× | **0.34 %** |

Distribution **très déséquilibrée et à longue traîne** (médiane ≪ moyenne ; quelques minimizers très
fréquents). **Mais** le plus gros bucket ne représente que **0.2–0.3 % du travail total** ⇒ avec un
**ordonnancement dynamique** (file de travail, ou buckets triés par taille décroissante), le parallélisme par
bucket atteint un makespan ≈ travail/P **quasi-linéaire** même à 22 threads. Un **découpage statique par plages
d'ID de bucket souffrirait du skew** (une plage peut concentrer plusieurs gros buckets).

---

## 7. Liste détaillée des points de ralentissement (classés par coût)

| # | Goulot | Coût (matérialisé) | Où | Nature |
|--:|---|--:|---|---|
| **1** | **Re-compaction de la sortie en super-k-mers** (`generate_sorted_list_from_enumeration` ré-exécuté par bucket) | **81–88 %** | `VirtualSkmer.hpp:387`, appelé `SetOperations.hpp:188` | algorithmique : on **jette l'ordre déjà produit par le merge** et on ré-assemble de zéro |
| **2** | **`colinear_chaining`** (tri + arbre de Fenwick + compression de coordonnées) par paire de colonnes par bucket | **~30 %** | `ColinearChaining.cpp` | travail + stall mémoire (1.18 cyc/inst) ; calcule une compaction **optimale** alors qu'une compaction valide suffit |
| **3** | **`get_candidate_overlaps`** (tri de `right_keys` + `lower_bound` par élément) | **~27 %** | `VirtualSkmer.hpp:627` | re-trie une colonne **déjà triée** par le merge |
| **4** | **`sort_column`** (`std::sort`+`std::unique`+pré-scan `has_valid_kmer`) par colonne par bucket | **~14 %** | `VirtualSkmer.hpp:585` | re-trie des **runs déjà ordonnés** ; dense mais rapide (0.53 cyc/inst) |
| **5** | **Merge driver + collecte** : `merge_columns` + `CollectSink::*` + `get_skmer_of_kmer` (copie par valeur) | **~10–15 %** | `SetOperations.hpp:56`, `Skmer.hpp:1109` | le plus stallé (1.58 cyc/inst) : sauts de curseur + matérialisation de skmers mono-k-mer |
| 6 | `merge_LList_column` (4e étape de l'assemblage) | ~6 % | `VirtualSkmer.hpp:678` | — |
| 7 | **Hole-scan du merge** : `set_next_valid` re-balaie tout depuis 0 à **chaque** colonne (k−m+1 passes) | ~31 % du **`_size`** (≈4 % du matér.) | `SetOperations.hpp:44,61` | redondance : chaque record est testé (k−m+1)× |
| 8 | **`kmer_compare`** copie le `pair` + masque 2× par appel | ~52 % du **`_size`** | `Skmer.hpp:1071` | micro : copie évitable |
| 9 | Largeur de record : **16 octets pour ~50 bits utiles** (+ padding) ; bascule uint128 à k>32 | bande passante du scan ; +15 % à k=41 | `Skmer.hpp:25` | layout |
| 10 | **Séquentiel** : un seul cœur ; buckets indépendants inexploités | facteur P perdu | `SetOperations.hpp:136,176` | parallélisme absent |
| 11 | I/O **cache chaud ≈ 0**, mais **cache froid** lit ~2× la taille des fichiers (chr1 : ~1.8 Go) → peut dominer le **`_size`** sur stockage lent | 0 (chaud) | `VirtualSkmer.hpp:1275` | non recouvert |
| 12 | Surcoût fixe par bucket (manipulateur + masques + `SortedVirtualSkmerList sub` recréés par bucket) | <1 % (gros buckets) ; ↑ si beaucoup de petits buckets | `SetOperations.hpp:171,188` | médiane ecoli = 6 records/bucket |

---

## 8. Liste exhaustive de suggestions d'accélération

Classées par retour sur investissement, avec impact attendu rattaché aux mesures.

### A. Algorithmique — attaquer les 85 % (re-compaction). Les gros gains.

- **A1 — Assembler les super-k-mers EN STREAMING pendant le merge, au lieu de collecter puis ré-assembler.**
  Le merge produit déjà les k-mers retenus **triés par colonne** ; au lieu de les jeter dans `collected` puis
  de relancer tri+overlaps+chaînage (goulots #2/#3/#4), étendre gloutonnement le super-k-mer courant tant que
  le k-mer suivant en est le successeur contigu. Supprime #2+#3+#4 d'un coup.
  **Impact : viser un temps matérialisé ≈ merge ⇒ ~3–5×** sur ∩/∪. *Le plus gros levier.*

- **A2 — Mode « sans compaction » (`--no-compact`).** Le résultat n'a besoin d'être qu'une liste triée valide
  et interrogeable. Écrire chaque k-mer résultat comme skmer mono-k-mer (sauter entièrement
  `generate_sorted_list_from_enumeration`) ⇒ matérialisé ≈ `_size` + écriture.
  **Impact : ~10×**, au prix d'un fichier ~4–5× plus gros (un record/k-mer). Idéal pour pipelines qui
  re-interrogent ou re-dérivent.

- **A3 — Exploiter l'ordre du merge dans la re-compaction (si A1 trop intrusif).**
  (i) `sort_column` (#4) opère sur des **runs déjà triés** → remplacer `std::sort` par un garde `is_sorted`
  ou une fusion k-voies de runs. (ii) `get_candidate_overlaps` (#3) re-trie `right_keys` **déjà triés** → le
  tri + `lower_bound` deviennent un balayage linéaire à deux curseurs. **Impact : la quasi-totalité de
  #4 (14 %) + le tri interne de #3.**

- **A4 — Remplacer le chaînage colinéaire optimal par un chaînage glouton O(n).** `colinear_chaining` calcule
  un **maximum** d'overlaps compatibles (tri + Fenwick + compression de coordonnées). Un chaînage glouton
  gauche→droite produit des super-k-mers **valides** (quasi-maximaux) sans tri ni arbre.
  **Impact : ~30 % du temps total** récupéré ; sortie marginalement moins compacte.

### B. Parallélisme — buckets indépendants (le meilleur ROI d'ingénierie après A)

- **B1 — Paralléliser la boucle par bucket** de `set_sizes`/`materialize_setop`. Les buckets sont totalement
  indépendants ; le plus gros = 0.2–0.3 % du travail (§6) ⇒ **ordonnancement dynamique** (file de travail, ou
  buckets triés par taille décroissante) pour un scaling **quasi-linéaire**. **Éviter** le découpage statique
  par plages d'ID (skew 10–14×). **Impact : ≈ T/P** (≈8× sur 8 cœurs) sur le chemin matérialisé dominant.
- **B2 — Matérialisation en deux passes pour une écriture parallèle sans verrou.** Passe 1 = merge `_size`
  parallèle ⇒ comptes par bucket ⇒ offsets fichier ; passe 2 = re-compaction+`pwrite` de chaque bucket à son
  offset. Supprime la sérialisation de l'écriture. (Combinable avec A1/A4.)
- **B3 — SIMD/vectorisation du merge** `_size` (cardinalité) : comparaisons de k-mers par lots.

### C. Micro-optimisations du merge (utile surtout pour `_size`/Jaccard, où le merge EST le coût)

- **C1 — `kmer_compare` sans copie** (`Skmer.hpp:1071`) : comparer via masques sur référence, et pour les
  largeurs uint32/uint64 comparer le mot directement plutôt que la struct `pair`. La copie + `&=` du `pair`
  pèsent ~15 % des instructions du `_size`.
- **C2 — Supprimer le re-balayage par colonne** dans `merge_columns` (#7). Indexer une fois chaque liste par
  colonne (CSR : pour chaque colonne, la liste de ses records) ⇒ le scan (k−m+1)×N devient N. Le m-sweep
  montre **jusqu'à −37 %** sur le merge. (Alternative : un seul passage entrelacé qui avance toutes les
  colonnes ensemble.)
- **C3 — Itérer seulement sur les colonnes valides** via `get_valid_kmer_bounds` (O(1) par record) au lieu de
  tester `has_valid_kmer` pour chaque (record, colonne). Élimine l'essentiel des 12 % de `has_valid_kmer`.

### D. Layout mémoire / bande passante

- **D1 — Compacter le record** (16 o pour ~50 bits utiles + 2×uint16 + padding) : replier les tailles
  pref/suff dans les bits libres du `pair`, ou utiliser la largeur exacte ⇒ **−25 % d'octets scannés** (aide
  le merge et le cache ; se cumule avec C). Surveiller le seuil uint128 à k>32.
- **D2 — Réutiliser les scratch entre buckets** : `collected` et les vecteurs internes de `sub` sont recréés
  par bucket (le scratch de fusion `m_merge_scratch` l'est déjà). Pré-dimensionner depuis la taille du bucket.
  Mineur (alloc < 1 %), mais gratuit.

### E. I/O (invisible à chaud ; réel à froid / stockage lent)

- **E1 — Recouvrir lecture et calcul** : `posix_fadvise(SEQUENTIAL|WILLNEED)` + readahead sur le flux de
  buckets, et **double-buffer** (lire le bucket i+1 pendant le calcul du bucket i). Sur le chemin `_size`
  (~3 s de calcul à chr1), une lecture à froid de ~1.8 Go peut dominer sur SSD SATA/HDD.
- **E2 — `mmap` du payload** pour laisser le noyau gérer le readahead ; bufferiser les grosses écritures
  (union de chr1 ≈ jusqu'à 380M k-mers).

### F. API / guidage d'usage (sans changement de code)

- **F1 — Documenter les variantes `_size`** : si seule la cardinalité (Jaccard/containment) est voulue, elles
  évitent **85–91 %** du travail. Déjà présentes ; à mettre en avant.
- **F2 — `--buckets 16384`** pour les charges set-op (≈19 % plus rapide que 4096 ; au-delà, rendements
  décroissants + surveiller la largeur à très grand nombre de buckets / grand k).
- **F3 — m libre** : choisir m pour la taille d'index / la requête, pas pour les set-ops (effet négligeable).

### Ce qui n'est **pas** un goulot (pour cadrer l'effort)

I/O à chaud (~0 %), allocation (<1 %), génération des masques / setup manipulateur (<1 % sur gros buckets),
largeur de record pour k≤31, et — malgré son **nombre** d'appels élevé — `has_valid_kmer` en **temps**
(cache-résident, ~quelques cycles). Le comptage d'appels surévalue le scan ; c'est `kmer_compare` et surtout
la **re-compaction** qui coûtent.

---

## 9. Reproductibilité

- Harnais : `scripts/bench/bottleneck/decompose.sh` (décomposition `_size` vs matérialisé) et
  `scripts/bench/bottleneck/perf_categorize.py` (agrégation `perf report` par phase).
- Données brutes versionnées : `report/data/bottleneck_decompose.csv`, `bottleneck_perf_summary.txt`,
  `bottleneck_sweep_{m,buckets,k}.csv`, `bottleneck_construct.csv`.
- Builds (hors dépôt, dans `/tmp`) : Release `-O3 -march=native -g -fno-omit-frame-pointer` pour le timing +
  `perf` ; `CMAKE_BUILD_TYPE=Profile` (`-fno-inline`) pour `callgrind`.
- Génomes : `…/sklib/scripts/out/e2e/genomes/*.sanitized.fa` (réels) ; copies mutées via
  `scripts/bench/mutate.py`.
- Exemples :
  - décomposition : `bash scripts/bench/bottleneck/decompose.sh`
  - perf : `perf record -F 2500 -g -o p.data -- sskm setop --op intersection -a A.sskm -b B.sskm -o /dev/null`
    puis `python3 scripts/bench/bottleneck/perf_categorize.py p.data`
  - callgrind (build Profile) : `valgrind --tool=callgrind --cache-sim=no sskm setop --op intersection_size -a A -b B`

---

## 10. Résultats d'implémentation (2026-06-07)

Les suggestions ont été implémentées **une par une** (sauf B1, exclu), chacune validée — **correction exacte**
(résultat vérifié contre un binaire *oracle* immuable : sous-ensembles + cardinalités, invariant au
re-packing ; + identité de requête contre vérité-terrain) puis **timing grande échelle** — et **gardée
seulement si le gain est réel**. Suite gtest complète après tous les changements : **196/197 OK** (1 skip
préexistant). Mesures cache chaud, 1 cœur, médianes.

### Gain cumulé (original → A2+A3+C2+D2+A4)

| opération | original | final | **speedup** |
|---|--:|--:|--:|
| celegans ∩ (compact) | 13.88 s | 8.74 s | **1.59×** |
| celegans ∪ (compact) | 20.32 s | 12.64 s | **1.61×** |
| chr1 ∩ (compact, 190M) | 35.01 s | 22.18 s | **1.58×** |
| celegans ∩ `_size` (cardinalité) | 2.06 s | 1.37 s | **1.51×** |
| chr1 ∩ `_size` (cardinalité, 190M) | 3.00 s | 2.21 s | **1.36×** |
| **celegans ∩ `--no-compact`** (nouveau mode) | 13.86 s | **4.63 s** | **2.99×** |

### Détail par amélioration (gardée ✅ / revert ⛔ / différée ⏸)

| # | amélioration | verdict | gain mesuré | correction |
|---|---|:--:|---|---|
| **A2** | mode `--no-compact` (écrit 1 record/k-mer, saute la re-compaction) | ✅ | **2.75–3.27×** (sortie ~3.3–5× plus grosse) | requête-identique vs vérité-terrain |
| **A4** | `greedy_chaining` (LIS *patience-sort*) au lieu du Fenwick `colinear_chaining`, set-op only | ✅ | **1.34–1.38×** (compaction identique) | même nb records ; ordre par colonne préservé (chaîne strictement croissante non-croisée) ; requête-identique |
| **C2** | index CSR par colonne dans `merge_columns` (supprime le rescan (k−m+1)×) | ✅ | **1.38–1.52× sur `_size`**, +5–7% matérialisé | **byte-identique** |
| **A3** | garde `is_sorted` dans `sort_column` (l'énumération set-op est pré-triée par colonne) | ✅ | +7–8% matérialisé | **byte-identique** (construction incluse) |
| **D2** | réutiliser le worker de re-compaction entre buckets (manipulateur/masques construits 1×) | ✅ | ~2% matérialisé | **byte-identique** |
| C1 | `kmer_compare` sans copie | ⛔ | 0.99× (aucun, `-O3` déjà optimal) | byte-identique |
| A1 | assemblage *streaming* pendant le merge | ⏸ | version contenue 0.98× ; complète différée | — |
| D1 | record compact (16→12 o) | ⏸ | non implémenté | — |

### Notes

- **A2** est le plus gros levier (×2.75–3.27) : à utiliser quand un index de sortie plus volumineux est
  acceptable (pipelines qui re-interrogent/re-dérivent). La sortie reste une liste triée valide et
  **s'interroge à l'identique** (vérifié contre vérité-terrain, pas seulement vs l'ancien code).
- **A4** : le `colinear_chaining` optimal (tri + arbre de Fenwick + compression de coordonnées) est
  remplacé, **sur le seul chemin set-op**, par une LIS *patience-sort* (un tri + une recherche binaire par
  overlap, sans Fenwick ni compression) — **chaîne maximale donc compaction identique**. La construction
  garde `colinear_chaining` (ses 30 tests restent verts). L'ordre par colonne est préservé car la chaîne
  est strictement croissante dans les deux coordonnées (propriété non-croisée, identique à colinear).
- **C2** rend la **cardinalité (`_size`, Jaccard) ~1.5× plus rapide**, byte-identique.
- **A1 différée** : après A3/A4/C2/D2 la re-compaction est déjà ~1.6× plus rapide ; le goulot résiduel
  (`get_candidate_overlaps`, ~35%) est un *join* sort-merge inhérent (borné cache, pas allocation) et
  **partagé avec la construction** (tests figés) → réécriture *streaming* à fort risque pour un gain
  marginal. La version contenue (réutilisation de buffers) ne gagne rien (alloc non-goulot).
- **D1 différée** : changerait le **format on-disk** (casse la compatibilité des `.sskm` existants) +
  grande surface (sérialisation/lecteurs/requête/déterminisme `m_pad`) pour ~quelques % (le merge ne pèse
  que 10–15% du matérialisé). Projet séparé à format versionné.
- **B1** (parallélisme par bucket, exclu à la demande) reste **le plus gros levier inexploité** :
  ordonnancement dynamique → ~T/P quasi-linéaire (cf. §6).

Reproductibilité de la validation : harnais `verify.sh` (correction, oracle) et `bench.sh` (timing) du
plan de travail ; toutes les améliorations gardées sont byte-identiques sauf A2/A4 (set-équivalentes,
requête-identiques, vérifiées contre vérité-terrain).
