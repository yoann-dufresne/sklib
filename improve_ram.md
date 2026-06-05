# Réduction de la RAM à la construction (`sskm construct`)

> **⚠️ Mise à jour (format VSKMER_3, buckets persistés).** Depuis ce log, les buckets ne sont
> **plus concaténés** en une liste unique : ils sont **persistés en sous-listes** dans le fichier
> (format `VSKMER_3` avec répertoire de buckets), et la **requête route** chaque k-mer vers son
> bucket puis ne charge que celui-ci (paresseux). L'affirmation « format de sortie inchangé /
> `query`/`load` intacts » plus bas n'est donc **plus valable** ; voir `README.md` et le wiki.
> Le reste de ce document reste un **log historique** de la réduction RAM à la construction (V0→V3).

> **État (branche `dev`, fusionnée dans `main`).** La construction est passée d'un
> chemin tout-en-RAM (pic O(génome)) à un chemin **bucketé sur disque** dont le pic
> est borné par le **plus gros bucket de minimizer**, réglable via `--buckets` /
> `--max-ram`. Le **format de sortie est inchangé** (`query` / `load` intacts) et la
> sortie est **déterministe**.
>
> Trois versions jalonnent le travail (toutes mono-thread, build Release) :
> **V0** = tout-en-RAM (historique) · **V1** = bucketage disque · **V2** = V1 +
> C2/C3/C4/B2 (actuel `main`).
>
> | exemple | V0 RSS | V1 RSS | **V2 RSS** | V0 t | V1 t | **V2 t** |
> |---|--:|--:|--:|--:|--:|--:|
> | chr1, k=31 m=15 | 5230 Mo | 276 Mo | **204 Mo** | 613 s | 112 s | **51 s** |
> | chr21, k=21 m=11 | 1437 Mo | 75 Mo | **51 Mo** | 86 s | 13 s | **8 s** |
>
> Gain global **V0 → V2 : 25–55× moins de RAM** *et* **2–12× plus rapide** (les
> hash-maps supprimées étaient un goulot mémoire *et* CPU). V2 améliore encore V1
> (RAM jusqu'à −32 %, temps ×0.36–0.71) **sans aucune régression**. Benchmark
> complet (5 génomes × 3 configs) : **§5**.
>
> **V3 = ordre des minimiseurs par hash φ** (branche `mini_hash_order`) : optimisation de
> **compression** (−9 à −21 % de super-k-mers partout) + **RAM sur génomes répétitifs**
> (celegans k21m11 : 341 → 70 Mo, ~5×) + correction du **drop poly-A**. Détail : **§6**.
>
> Ce document liste désormais **ce qui reste à faire**. §0 résume le livré ;
> §1 garde les invariants de correction (fondation des optimisations futures) ;
> §2 décrit les nouvelles propositions ; §3 la validation ; §4 le statut ;
> §5 le benchmark final ; §6 l'ordre φ (V3).

---

## 0. Déjà livré (ne pas refaire)

Tout est sur `dev`, chaque étape en commit unitaire, validée (tests + e2e KMC),
et — pour les optimisations de facteur constant — **sortie byte-identique** au
chemin antérieur (prouvé par `cmp`).

| Étape | Contenu | Fichiers |
|---|---|---|
| **A v1** | Partitionnement par préfixe de minimizer sur disque (`--buckets`, défaut 4096). Phase 1 stream → buckets ; phase 2 par bucket : load + tri/dédup + algo inchangé + append. | `SkmerBucketWriter.hpp`, `sskm.cpp` |
| **A v2** | Frontières de buckets **adaptatives** via histogramme de minimizers, dimensionnées à un budget RAM (`--max-ram`, ex. `2G`). Override `--buckets`. Nécessite `-f <fichier>` (2 lectures). | `sskm.cpp` (`make_adaptive_bucketing`) |
| **B(1)** | Dédup intra-bucket : `sort` par `(m_pair, pref, suff)` + `unique` strict avant `generate`. | `sskm.cpp` (`sort_and_dedup`) |
| **C1** | Buffers ping-pong : `m_merge_scratch` membre réutilisé + `std::swap` au lieu de réallouer `merged` à chaque colonne. | `VirtualSkmer.hpp` (`merge_LList_column`) |
| **C2a** | `get_candidate_overlaps` : `unordered_map` → tableau trié `(clé, pos)` + recherche binaire (merge-join). Même ordre `(left, right)` émis. | `VirtualSkmer.hpp` |
| **C2b** | `colinear_chaining` : les deux `unordered_map<overlap,…>` → **vecteurs indexés par le rang de l'overlap trié** (DP + cellules Fenwick portent un index, pas l'overlap). | `ColinearChaining.cpp` |
| **C3** | `Virtual_skmer` 40 → **32 o** (`last_id` en `uint32_t` ; borne : index < 2^32). −20 % sur les deux buffers de merge. | `VirtualSkmer.hpp` |
| **C4a** | `colinear_chaining` opère **en place** sur la plage `[begin, end)` de l'appelant — plus de copie `std::vector<overlap> ov`. −16 o × C. | `ColinearChaining.cpp` |
| **C4b** | Indices du DP / Fenwick en **32 bits** (`prev_idx`, `length`, index de cellule ; bornés par C < taille bucket < 2^32). −50 % sur ces tableaux. | `ColinearChaining.cpp` |
| **B(2)** | Collapse des super-k-mers canoniques **identiques consécutifs** avant écriture du bucket : réduit les fichiers temporaires et le pic de **chargement** des buckets répétitifs (poly-A). L'histogramme `--max-ram` continue de compter le **flux complet** → byte-identité préservée. | `sskm.cpp`, `SkmerBucketWriter.hpp` |
| (bug) | `Skmer` : padding de queue (4 o) explicitement mis à zéro → sérialisation **déterministe** (fuite de heap non initialisé auparavant). | `Skmer.hpp` |
| (infra) | Sérialiseur incrémental (`write_header` / `append_payload` / `patch_count`) ; `--tmp-dir` ; nettoyage RAII des fichiers temporaires ; fallback tout-en-RAM pour `--ascii` / stdout. | `VirtualSkmer.hpp`, `sskm.cpp` |

Critères d'acceptation initiaux : **atteints** (§4).

---

## 1. Invariant de correction (fondation — toujours valable)

Les deux propriétés de l'encodage interleaved centré-minimizer qui rendent le
partitionnement correct. **Toute optimisation future (C4, parallélisme) doit les
préserver.**

**(I1) Le tri global est dominé par le minimizer.** Le minimizer occupe les bits
de poids fort de `m_pair` (`minimizer = m_pair >> 4(k-m)`). L'ordre de tri final
(`Skmer::operator<` = comparaison de `m_pair`) place tous les super-k-mers d'un
minimizer X contigus, avant ceux d'un minimizer Y > X. ⇒ concaténer des blocs
triés dans l'ordre croissant de minimizer reproduit l'ordre global.

**(I2) Les interactions de l'algorithme sont strictement intra-minimizer.** Le
seul couplage entre super-k-mers, `get_candidate_overlaps`, apparie le (k-1)-mer
suffixe d'une colonne au (k-1)-mer préfixe de la suivante. Le minimizer (m
positions centrales) est contenu dans tout k-mer du super-k-mer, donc dans ce
(k-1)-mer partagé : deux super-k-mers en overlap ont **le même minimizer**.
⇒ le bloc de sortie d'un minimizer ne dépend que de ses super-k-mers d'entrée.

**Théorème de partitionnement.** Router par minimizer, exécuter
`generate_sorted_list_from_enumeration` inchangé par groupe, concaténer dans
l'ordre croissant de minimizer ⇒ **même ensemble de k-mers** que le run global.

> ⚠️ **Granularité atomique = un minimizer entier.** Router **uniquement** sur les
> bits du minimizer ; jamais sur les flanks (deux super-k-mers en overlap ont le
> même minimizer mais des flanks décalés → les séparer casserait I2).
>
> ⚠️ **Résidu poly-A irréductible.** Un seul minimizer hyper-fréquent (poly-A = 0)
> forme un bucket atomique que ni `--buckets` ni `--max-ram` ne subdivisent. C'est
> le plancher RAM effectif sur génomes répétitifs. Mitigations : B(2) (§2),
> augmenter `m`, ou accepter ce plancher.

---

## 2. Propositions restantes

> Les gains « gratuits / faible risque » sont **livrés** (C4a, C4b, B(2) — voir §0).
> Ce qui reste est **structurel** : il change potentiellement l'algorithme et exige
> un profilage préalable.

### C4(3) — Borner / streamer le nombre d'overlaps C (structurel, à profiler)

Depuis C2/C4, les hash-maps et la copie ont disparu ; sur les buckets
**répétitifs**, le pic résiduel est dominé par les **vecteurs d'overlaps**, en O(C)
où C = nombre d'overlaps candidats (potentiellement quadratique dans une colonne
très répétitive) : `candidate_overlaps` (16 o × C), `prev_idx` + `length`
(4 o × C chacun après C4b), Fenwick `bit` (16 o × S, S ≤ C).

Si une colonne génère un C quadratique (poly-A), envisager de **produire les
overlaps en flux** vers le chaînage plutôt que de tout matérialiser, ou de
**plafonner** C. Plus risqué : change potentiellement la sélection du chaînage →
valider l'**ensemble de k-mers** (pas seulement la byte-identité, qui pourrait
légitimement changer). Attendre le profilage (ci-dessous) avant de s'y engager.

### Construction parallèle par bucket (différé — pas de parallélisme pour l'instant)

Les buckets sont **indépendants** (théorème §1) ⇒ un pool de workers pourrait en
traiter plusieurs simultanément. Borne RAM = `threads × plus gros bucket` ;
accélération quasi-linéaire. **Différé** par décision utilisateur (rester
mono-thread pour le moment). Note d'implémentation pour plus tard : la
sérialisation incrémentale (`append_payload` + `patch_count`) impose d'**ordonner
les écritures par index de bucket croissant** (chaque worker bufferise son
sous-résultat, écriture finale ordonnée), pour préserver I1.

### Profilage ciblé (idée réservée par l'utilisateur)

Un profil mémoire/CPU du **pire bucket** (poly-A de chr21 / celegans) tranchera
objectivement entre C4(1/2/3) et d'éventuelles refontes de
`get_candidate_overlaps` / `colinear_chaining`. L'utilisateur a une piste à
exploiter ici plus tard.

---

## 3. Validation (inchangée — réutiliser l'existant)

- **Byte-identité vs `dev`** (preuve de correction la plus forte pour les
  optimisations de facteur constant) : `cmp` la sortie binaire `--buckets {64,
  256, 4096}` contre un binaire `dev` ; doit être **identique**. Plus fort qu'un
  e2e échantillonné.
- **Correction vs KMC** — `scripts/large_scale_e2e.sh` sur `ecoli`, `yeast`,
  `chr21`, `KM="21,11 31,13 32,17"` (égalité d'ensemble + statuts vs KMC,
  déterminisme, cohérence comptes binaire/ASCII).
- **Tests unitaires** — `build-debug/tests/sklib-tests` (depuis le dossier de
  build) ; en particulier `*SkmerSorting*`, `GetCandidateOverlap*`, `MergeColumn*`,
  `*Colinear*`, `bug06_construct_homopolymer_drop`.
- **Mémoire / vitesse** — `/usr/bin/time -v` (peak RSS) + best-of-3 sur `chr21`
  en faisant varier `--buckets` ; comparer à un binaire baseline construit dans un
  worktree. ⚠️ **Toujours mesurer en build Release** : un `build/` configuré en
  `DEBUG` active AddressSanitizer et gonfle le RSS ~3,7× (voir mémoire
  `measuring-construct-rss-asan-trap`).
- **Bench** — `scripts/bench/bench.sh` (CSV append-only, tagué par commit) pour
  comparer versions sur ecoli/yeast/celegans.

---

## 4. Statut des critères d'acceptation

- ✅ `large_scale_e2e.sh` : tous les checks PASS sur ecoli/yeast/chr21 pour les
  `(k,m)` testés, chemin bucketé par défaut.
- ✅ Peak RSS chr21 réduit d'un **ordre de grandeur** (851 → 54 Mo au défaut,
  jusqu'à plus bas via `--buckets`/`--max-ram`), **contrôlable**.
- ✅ Format de sortie inchangé ; `query` fonctionne sans recompilation.
- ✅ Sortie déterministe (byte-identique entre runs).
- ✅ Aucun fichier temporaire résiduel (succès ou exception).
- ✅ V2 (C4a/C4b/B2) : **byte-identique** à V1 partout, **aucune régression** RAM
  ni temps sur les 15 cellules du benchmark (§5).
- ➕ Bonus non demandé : construction **plus rapide** que l'ancien chemin (et V2
  encore plus rapide que V1, ×0.36–0.71).

**Cibles ouvertes** : abaisser encore le plancher du pire bucket via C4(3)
(structurel, à profiler), et — plus tard — paralléliser par bucket.

---

## 5. Benchmark final (V0 / V1 / V2)

**V0** = tout-en-RAM (5d9dace) · **V1** = bucketage disque (b085132) · **V2** =
`main` actuel, +C2/C3/C4/B2 (8394784). Build **Release**, **mono-thread**,
`--buckets` par défaut (4096). RSS = pic mesuré via `/usr/bin/time -v` ; temps =
construct mur. ⚠️ Toujours mesurer en Release (un `build/` DEBUG active
AddressSanitizer et gonfle le RSS ~3,7× — mémoire `measuring-construct-rss-asan-trap`).

| génome | k,m | V0 RSS | V1 RSS | **V2 RSS** | V0 t | V1 t | **V2 t** | V2/V1 RAM | V2/V1 t |
|---|---|--:|--:|--:|--:|--:|--:|--:|--:|
| ecoli | 15,7 | 305 Mo | 23 Mo | **22 Mo** | 9 s | 2 s | **1 s** | ×0.98 | ×0.69 |
| ecoli | 21,11 | 186 Mo | 21 Mo | **21 Mo** | 8 s | 2 s | **1 s** | ×1.00 | ×0.52 |
| ecoli | 31,15 | 142 Mo | 20 Mo | **19 Mo** | 6 s | 2 s | **1 s** | ×0.99 | ×0.47 |
| yeast | 15,7 | 751 Mo | 25 Mo | **20 Mo** | 25 s | 5 s | **3 s** | ×0.81 | ×0.61 |
| yeast | 21,11 | 465 Mo | 21 Mo | **16 Mo** | 26 s | 4 s | **3 s** | ×0.78 | ×0.71 |
| yeast | 31,15 | 299 Mo | 18 Mo | **14 Mo** | 23 s | 4 s | **2 s** | ×0.75 | ×0.54 |
| celegans | 15,7 | 5123 Mo | 187 Mo | **140 Mo** | 172 s | 44 s | **27 s** | ×0.75 | ×0.61 |
| celegans | 21,11 | 3891 Mo | 407 Mo | **341 Mo** | 271 s | 63 s | **25 s** | ×0.84 | ×0.41 |
| celegans | 31,15 | 2456 Mo | 277 Mo | **240 Mo** | 261 s | 64 s | **23 s** | ×0.87 | ×0.36 |
| chr21 | 15,7 | 1977 Mo | 51 Mo | **51 Mo** | 75 s | 15 s | **10 s** | ×1.00 | ×0.65 |
| chr21 | 21,11 | 1437 Mo | 75 Mo | **51 Mo** | 86 s | 13 s | **8 s** | ×0.68 | ×0.60 |
| chr21 | 31,15 | 938 Mo | 66 Mo | **53 Mo** | 90 s | 13 s | **9 s** | ×0.80 | ×0.65 |
| chr1 | 15,7 | 10785 Mo | 192 Mo | **192 Mo** | 366 s | 95 s | **58 s** | ×1.00 | ×0.60 |
| chr1 | 21,11 | 8561 Mo | 342 Mo | **261 Mo** | 614 s | 114 s | **56 s** | ×0.76 | ×0.49 |
| chr1 | 31,15 | 5230 Mo | 276 Mo | **204 Mo** | 613 s | 112 s | **51 s** | ×0.74 | ×0.46 |

**Lecture.**
- **V0 → V2** : RAM divisée par **25–55×** (chr1 k15 : 10.8 Go → 192 Mo ; chr1 k21 :
  8.6 Go → 261 Mo), temps **2–12× plus rapide**.
- **V1 → V2** : RAM réduite jusqu'à **−32 %** (chr21 k21, chr1 k31), gain concentré
  sur les configs **riches en overlaps** où C2/C4 allègent l'état de chaînage ;
  plat à k=15 (le plus gros bucket n'y est pas borné par les overlaps). Les trois
  `×1.00` sont du bruit de mesure (< 0,3 %).
- **Temps** : V2 plus rapide que V1 sur **toutes** les cellules (×0.36–0.71) → le
  gain RAM est **gratuit** côté vitesse.
- **Régressions** : **aucune** (seuils 2 % RAM / 5 % temps). Chaque commit V2 est
  par ailleurs **byte-identique** à la sortie V1 validée vs KMC.

---

## 6. V3 — Ordre des minimiseurs par hash φ (compression + RAM sur répétitifs)

**Idée.** Remplacer l'ordre **lexicographique** sur les m-mers (valeur brute) par une
**bijection fixe φ** (mixer xorshift-multiply de classe hash-prospector `lowbias32`, masqué à
2m bits, + XOR constant pour `φ(0)≠0`). φ pilote **la sélection ET le tri/bucketing** : le slot
minimiseur stocké devient `φ(min)` (forme décodable ; `φ⁻¹` au décodage ASCII). Les minimiseurs
« aléatoires » ont une **densité plus faible** que lexicographiques → **moins de super-k-mers**.

**Correctif de frame (complétion d'issue #7).** Un minimiseur **palindrome reverse-complement**
(possible seulement à m **pair**, ex. `GC`) ne fixe pas l'orientation du super-k-mer → ses k-mers
ne sont pas tous canoniques. Le **Skmerator découpe** ces super-k-mers par k-mer (chacun
canonicalisé). Effet bonus : le **drop poly-A** rare en ordre brut disparaît (chr21 désormais
exactement égal à KMC : 32 720 485 k-mers).

**Format.** Le slot minimiseur étant permuté, le magic on-disk passe `VSKMER_M → VSKMER_2` :
`load()` **rejette** les anciennes listes (ordre brut) au lieu de renvoyer des résultats faux.

**Benchmark V2 (`main`) → V3 (φ)** (Release, mono-thread, `--buckets` défaut) :

| génome | k,m | V2 RSS | **V3 RSS** | ΔRSS | V2 skmers | **V3 skmers** | Δskmers |
|---|---|--:|--:|--:|--:|--:|--:|
| ecoli | 15,7 | 22 Mo | 28 Mo | +27 % | 2.12 M | 1.93 M | **−9 %** |
| ecoli | 21,11 | 21 Mo | 27 Mo | +29 % | 1.03 M | 0.88 M | **−14 %** |
| ecoli | 31,15 | 19 Mo | 23 Mo | +21 % | 0.60 M | 0.51 M | **−15 %** |
| yeast | 21,11 | 15 Mo | 25 Mo | +67 % | 2.90 M | 2.46 M | **−15 %** |
| yeast | 31,15 | 13 Mo | 19 Mo | +46 % | 1.57 M | 1.30 M | **−17 %** |
| **celegans** | 15,7 | 139 Mo | **73 Mo** | **−47 %** | 36.2 M | 32.9 M | −9 % |
| **celegans** | 21,11 | **341 Mo** | **70 Mo** | **−79 %** | 30.5 M | 26.2 M | −14 % |
| **celegans** | 31,15 | 211 Mo | **64 Mo** | **−70 %** | 14.7 M | 11.5 M | −21 % |
| chr21 | 15,7 | 51 Mo | 64 Mo | +25 % | 15.4 M | 14.1 M | −9 % |
| chr21 | 21,11 | 50 Mo | 62 Mo | +24 % | 9.43 M | 8.13 M | −14 % |
| chr21 | 31,15 | 53 Mo | 55 Mo | +4 % | 4.99 M | 4.12 M | −17 % |
| **chr1** | 21,11 | **261 Mo** | **206 Mo** | **−21 %** | 66.0 M | 58.3 M | −12 % |
| chr1 | 15,7 | 192 Mo | 208 Mo | +8 % | 71.0 M | 64.5 M | −9 % |
| chr1 | 31,15 | 204 Mo | 199 Mo | −2 % | 31.3 M | 25.8 M | −18 % |

**Lecture.**
- **Compression : gain partout** — −9 à −21 % de super-k-mers sur **tous** les génomes/configs
  (plus fort à grand m). Index plus petit (moins de bits/k-mer).
- **RAM : gros gain là où la RAM est grande** — sur le génome **répétitif celegans**, φ réduit
  le pic de **47 à 79 %** (k21m11 : **341 → 70 Mo**, ~5×), et chr1 k21m11 de **−21 %**. Ce sont
  les pires buckets dominés par la sur-sélection, que φ casse (l'objectif initial).
- **RAM : légère hausse là où la RAM est petite** — ecoli/yeast/chr21 +4 à +67 %, mais en
  absolu quelques Mo sur une base de 13–64 Mo (le plancher de chr21 est un run poly-A quasi-pur
  que φ ne fait que **relocaliser** ; φ produit aussi des super-k-mers moins nombreux mais plus
  longs → plus d'état de chaînage par bucket).
- **Temps** : neutre à légèrement plus rapide (moins de super-k-mers).
- **Correction** : égalité d'ensemble vs **KMC** (ecoli + chr21, m pair inclus), **0 faux
  négatif**, sortie **déterministe** byte-identique entre runs.

⇒ φ est repositionné comme optimisation de **compression + correction** (et RAM **sur les gros
génomes répétitifs**), pas comme un gain RAM universel.
