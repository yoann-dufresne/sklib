# OPTIM_SWEEP — campagne d'optimisation repartant de zéro (2026-07)

Exploration complète des performances de sklib **sans aucun prior** : aucun journal antérieur lu,
aucun historique git consulté, aucun chiffre repris. Tout ce qui suit vient du code source, du
harnais de benchmark et de mesures faites dans cette session.

- **Base de départ** : `572dc15` (v0.13.2, branche `dev`)
- **Branche livrée** : `optim/2026-07`, 7 commits
- **Machine** : Intel Core Ultra 7 165H — 6 P-cores SMT (cpu 0–11, 4,7–5,0 GHz), 8 E-cores
  (12–19, 3,8 GHz), 2 LP-E (20–21, 2,5 GHz) ; 62 GB RAM ; gouverneur `powersave` ; NVMe.
- **Compilateur** : clang++-18, `Release` = `-O3 -march=native` + LTO (déjà actif avant campagne).

## Résumé

| Workload | Régime | Gain |
|---|---|---|
| `construct` | k=95 (`kuint256`) | **+49 %** (t=1) · **+47 %** (t=6) |
| `construct` | k=15/31/63 | **+7 % à +11 %** |
| `query` | k=31, t=1 | **+37 % à +54 %** selon le taux de présence |
| `query` | k=63 / k=95, t=1 | **+24 à +28 %** / **+51 %** |
| `query` | k=31, multithread | **+172 %** (t=2) · **+23 %** (t=4) · **+18 %** (t=6) |
| `setop` | matérialisé, k=31/95 | **+2 à +3,5 %** |
| `setop` | `--sizes`, k=15/63 | neutre |

Sortie **byte-identical** à la baseline dans tous les cas : aucun changement de format, aucun
changement de valeur par défaut de la CLI, aucune perte d'exactitude.

---

## 1. Protocole et plancher de bruit

**C'est la partie qui a le plus changé mes conclusions.** Première tentative : séries A puis B,
consécutives, 9 répétitions chacune. Résultat sur `query chr21 k31 -t1` : médiane A = 1,330 s,
médiane B = 1,463 s — **10 % d'écart entre deux séries du même binaire exécutant la même
commande**, alors que l'IQR intra-série valait 2–3 %. La dérive thermique/turbo entre séries
dépassait largement les gains cherchés.

Protocole retenu, appliqué à **toutes** les mesures de ce rapport :

- **alternance A/B run par run** dans une seule boucle (`measure_ab`), ce qui annule la dérive ;
- pin sur des **P-cores physiques distincts** — `taskset -c 6` (t=1), `0,1,3,6,8,10` (t=6).
  Jamais `0-7`, qui apparie des siblings HT (cpu1/2, 3/4, 6/7… partagent un cœur) ;
- 5 à 9 répétitions, **médiane + IQR**, un run de chauffe jeté ;
- chaud par défaut (page cache amorcé) ; froid via `posix_fadvise(DONTNEED)`, noté quand utilisé.

**Plancher de bruit mesuré** (même binaire des deux côtés, protocole alterné) :

| Régime | B/A du témoin | IQR intra-série | Seuil de décision retenu |
|---|---|---|---|
| construct t=1 | 1,0052 | 0,95–1,19 % | **> 1,5 %** |
| construct t=6 | 1,0011 | 2,6–6,5 % | **> 5 %** |
| query t=1 | 1,0023 | 2,2–3,2 % | **> 1,5 %** |
| query t=6 | 1,0182 | 9,2–12,0 % | **> 5 %**, souvent non concluant |

Conséquence honnête : **à t=6 sur des charges courtes (< 1 s), rien en dessous de 5 % n'est
mesurable sur cette machine.** Plusieurs verdicts multithread ci-dessous sont donc « non
concluant » plutôt que « nul ».

Corpus figé : indexes construits une fois par (dataset, k, m) — ecoli 4,6 Mb, yeast 12 Mb,
chr21 39 Mb, chr1 224 Mb ; requêtes = 20 000 reads × 300 bp par taux de présence (0/50/100 %),
plus un jeu de 120 000 reads ; set-ops contre des mutants aux Jaccard mesurés 0,103 / 0,496 /
0,896. Échelle des largeurs de stockage (à `--buckets 4096`, donc b=12) :

| (k,m) | bits retenus `2(2k−m)−b` | store | octets/enregistrement |
|---|---|---|---|
| 15,7 | 34 | `uint32_t` | 16 (12 utiles) |
| 31,15 | 82 | `uint64_t` | 24 (20 utiles) |
| 63,31 | 178 | `__uint128_t` | 48 (36 utiles) |
| 95,47 | 274 | `kuint256` | 72 (68 utiles) |

---

## 2. Carte du temps (phase 2 — avant toute modification)

### construct — répartition phase 1 / phase 2 (chr21, `SKLIB_TIMING`)

| k | t=1 phase 1 | t=1 phase 2 | t=6 phase 1 | t=6 phase 2 |
|---|---|---|---|---|
| 15 | 33,0 % | 67,0 % | 28,5 % | 71,5 % |
| 31 | 38,9 % | 61,1 % | 35,9 % | 64,1 % |
| 63 | 46,3 % | 53,7 % | 41,8 % | 58,2 % |
| 95 | **74,8 %** | 25,2 % | **70,2 %** | 29,8 % |

Profil plat (`perf`, t=1, chr21) :

| k=31 | | k=63 | | k=95 | |
|---|---|---|---|---|---|
| `generate_sorted_list_from_enumeration` | 23,6 % | `SeqSkmerator::operator++` | 25,6 % | `minimizer_is_ambiguous` | 22,3 % |
| `SeqSkmerator::operator++` | 19,4 % | `generate_sorted…` | 21,3 % | `operator++` | 21,2 % |
| **`sort_column`** | **11,1 %** | `minimizer_is_ambiguous` | 11,9 % | **`permute_minimizer_slot`** | **15,8 %** |
| `greedy_chaining` | 7,6 % | **`sort_column`** | **11,7 %** | **`minimizer_rank`** | **12,7 %** |
| introsort de `sort_column` | 6,3 % | introsort de `sort_column` | 4,5 % | `generate_sorted…` | 10,5 % |
| introsort de `sort_and_dedup` | 3,8 % | `greedy_chaining` | 5,3 % | `sort_column` | 5,6 % |

### query — stream p50, chr1, t=1

| | k=31 | k=63 |
|---|---|---|
| `search_kmers_in_span_into` | **63,0 %** | **46,1 %** |
| producteur (`operator++` + `minimizer_is_ambiguous` + `init_seq`) | 11,3 % | 20,7 % |
| `__memset_avx2` (zero-fill de `resize` au chargement de bucket) | 4,1 % | 7,1 % |

Compteurs matériels (chr1 k=31, t=1) : IPC 1,57 ; 108 M cache-misses ; 10,6 M LLC-load-misses ;
**55,4 M branch-misses pour 5,4 M k-mers** (≈ 10 par k-mer — la dichotomie est imprévisible par
construction).

Attribution **par ligne source** (build Release + `-g`, même codegen) :

| ligne | % | quoi |
|---|---|---|
| `Skmer.hpp:1397` | **16,0 %** | `get_valid_kmer_bounds` — en réalité **le défaut de cache sur l'enregistrement sondé** |
| `VirtualSkmer.hpp:168-171` | 7,4 % | `update_searchable_positions` (balayage O(colonnes) par sonde) |
| `Skmer.hpp:1390` | 4,4 % | `has_valid_kmer` |
| `VirtualSkmer.hpp:349` | 3,7 % | site d'appel de `get_valid_kmer_bounds` |

### setop — union matérialisée, chr21 k=31 J=0,5, t=1

Décomposition wall (`SKLIB_UNION_PHASE_TIMING`) : read 1,4 % · **merge+collect 24,5 %** ·
**recompaction 71,8 %** · write 2,3 %.

`perf` (avec inlining) : `generate_sorted_list_from_enumeration` 52,2 % self / 70,8 % children ·
`greedy_chaining` 15,2 % · `materialize_setop` (merge + sink) 20,0 % · `build_column_csr` 4,5 %.

`--sizes` (compte seul) : `merge_columns` 53,6 % · `build_column_csr` 29,5 % ·
`mark_identical_records` 10,6 %.

Attribution par ligne dans la recompaction : **~9,6 % en réallocation de `std::vector`**
(`_M_realloc_insert`, `push_back`) — du churn pur, pas du travail de chaînage ; ~7 % dans la
jointure de hachage des overlaps ; ~5,8 % dans le prédicat de la boucle de merge.

---

## 3. Inventaire des leviers

19 leviers instruits. **7 retenus et livrés**, **7 testés et rejetés sur mesure**, **5 non testés**.

### Retenus (livrés)

| # | Levier | Commit | Gain mesuré |
|---|---|---|---|
| L1 | CSR de colonnes dans la compaction (construct) | `ab9baa3` | construct +5,9 % (k31) / +6,3 % (k63) |
| L12–L14 | Chemins rapides 128 bits pour 64 < 2m ≤ 128 | `508d3cf` | construct k95 **+60,9 %**, query k95 +49,3 % |
| L7 | Court-circuit XOR dès 64 bits (seuil 16 → 8 o) | `e30d64b` | query k31 +2,6 à +3,8 % |
| L18 | Réutilisation des buffers d'overlaps/chaîne | `a77e22c` | setop +3,3 % (union k31), +3,6 % (combiné) |
| L19 | Table de préfixe de minimizer par bucket | `76f6f17` | query **+22 à +55 %** |
| L4 | Parallélisation du *parse* de requête | `41ad75b` | query t=2 +173 %, t=4 +67 %, t=6 +44 % |
| L8 | Tri de colonne sur clés pré-extraites | `d3ec019` | construct +4,9 % (k15) / +2,8 % (k63) |

### Rejetés sur mesure

| # | Levier | Verdict chiffré | Pourquoi |
|---|---|---|---|
| L15 | Supprimer le zero-fill de `resize()` au chargement de bucket | +0,53 / +0,09 / +0,74 / −0,50 % (chaud, froid, t=1, t=6) | Les 4–7 % attribués à `__memset_avx2` par `perf` sont du **défaut de page déplacé**, pas du travail supprimé : les pages fraîchement allouées doivent être fautées de toute façon, la lecture les touche ensuite. Sous le plancher partout. |
| L17 | `get_skmer_of_kmer` par référence au lieu de par valeur | −0,33 % / −0,92 % | Le compilateur élidait déjà la copie (fonction inlinée, argument non modifié). |
| L16 | Mémoïser `col_skmer` dans `merge_LList_column` | +1,25 / −0,17 / +0,83 % | La redondance existe (le côté liste gagne ⇒ `col_idx` n'avance pas ⇒ recalcul), mais elle est trop rare pour dépasser le plancher. |
| L2 | **Packing des enregistrements (bump `VSKMER_6`)** | **inapplicable en l'état** | `__attribute__((packed))` est refusé par clang : `-Wpacked-non-pod`, *"not packing field 'm_pair' as it is non-POD for the purposes of layout"* — `Skmer::pair` a des constructeurs utilisateur. Détail en §5. |
| L10 | PGO (clang `-fprofile-generate/-use`) | **surapprentissage net** | +4,1 à +4,5 % sur les configurations vues à l'entraînement, mais **−6,6 % (setop k=63)** et **−9,8 % (construct k=63)** sur des configurations non vues. Le code étant massivement templaté par largeur, PGO optimise les instanciations profilées et pessimise les autres. |
| — | Conditionner le prefetch de la dichotomie à la largeur de plage | +1,01 / −0,34 / +2,05 % | Après L19 les plages sont minuscules et le prefetch coûte 4,5 % du profil, mais le gain mesuré est incohérent et sous le plancher. |
| L6 | E/S temporaires de phase 1 | **quantifié, non rentable** | `strace -c` : 52 635 syscalls, **0,175 s = 5,4 % du wall construct** (18 925 `openat`, 14 803 `write`, 18 907 `close`) — `flush_bucket` rouvre le fichier à chaque vidage. Mais : un cache de fd est inutile (4096 buckets touchés dans un ordre quasi aléatoire ⇒ taux de succès ~6 %), et doubler le budget d'écriture doublerait 32 Mo × n_threads de tampons. Au mieux ~3 % récupérables au prix de RAM. Détail en §5. |

### Non testés (et pourquoi)

| Axe | Candidat | Raison |
|---|---|---|
| Disposition | Cache de bucket en SoA (`hi[]`/`lo[]`/`sizes[]`) | Rendu largement caduc par L19 : la recherche ne touche plus que quelques enregistrements par requête. Le gain résiduel ne justifiait plus le coût. |
| E/S | `mmap` des payloads de bucket | Bloqué par l'alignement : `payload_start = 56 + 16·n_buckets ≡ 8 (mod 16)`, donc un pointeur zéro-copie serait mal aligné pour les enregistrements de 16/32 o. Débloquerait avec un bump de format (voir §5). |
| Mémoire | Huge pages (`MADV_HUGEPAGE`) sur les payloads résidents | `dTLB-load-misses` = 6,8 M contre 5,9 G cycles — sous le seuil d'intérêt. |
| Algorithmique | Suppression du balayage O(colonnes) par sonde (`update_searchable_positions`, 7,4 %) | L19 a réduit le nombre de sondes d'un facteur ~4, donc la part de ce balayage a fondu. À rouvrir si la requête redevient dominante. |
| Algorithmique | Recompaction des set-ops (toujours 71,8 % du coût matérialisé) | Le cœur (`merge_LList_column` + jointure d'overlaps + chaînage) est un algorithme, pas du gaspillage : les deux redondances identifiées (L16, L18) ont été traitées, L18 retenu. Une refonte demanderait de repenser la reconstruction des super-k-mers — hors périmètre d'une passe d'optimisation. |

---

## 4. Détail des leviers livrés

### L1 — CSR de colonnes dans la compaction (`ab9baa3`)

**Mécanisme.** `generate_sorted_list_from_enumeration` demandait à `sort_column` les ids valides
pour chacune des k−m+1 colonnes, et **chaque appel re-balayait toute l'énumération du bucket** en
testant `has_valid_kmer` sur chaque enregistrement : O(n × ncols) touches d'enregistrement, soit
n×ncols×sizeof(Skmer) octets de trafic (48 o × 33 colonnes à k=63). Or un enregistrement porte un
k-mer exactement sur l'intervalle contigu que renvoie `get_valid_kmer_bounds`, donc **une seule
passe** détermine l'appartenance de toutes les colonnes. Le CSR est construit par tableau de
différences (+1 à `start`, −1 après `end`, puis somme préfixe : O(1) par enregistrement au
comptage) puis une passe de placement. ~180 o/enregistrement au lieu de ~1584 à k=63.

**Identité.** La tranche contient exactement les ids que le balayage collectait, dans le même
ordre croissant ; les étapes `is_sorted`/`sort`/`unique` sont inchangées. `sort_column` reste en
place comme implémentation de référence épinglée par les tests `SortingColumn*`.

| Mesure (chr21, N=7) | Gain |
|---|---|
| construct k=15 t=1 | +2,38 % |
| construct k=31 t=1 / t=6 | +5,92 % / +5,29 % |
| construct k=63 t=1 / t=6 | +6,28 % / +6,23 % |
| construct chr1 k=31 t=6 | +4,33 % |

Croisé : query +0,49 %, setop −1,19 % (bruit ; le chemin set-op passe par `col_offsets` et saute
le CSR). Pic RSS 289,1 → 288,7 MB.

### L12–L14 — chemins rapides 128 bits pour 64 < 2m ≤ 128 (`508d3cf`)

**Mécanisme.** Trois routines par base / par yield avaient déjà un chemin « faire l'arithmétique
en `uint64_t` quand la valeur ne demande que 2m ≤ 64 bits, même si `kuint` est plus large ».
Au-delà de 2m = 64 elles retombaient toutes en `_BitInt(256)` — c'est-à-dire **toute la gamme
k ≈ 65…127**. `perf` à k=95 mettait 50,8 % du wall dans ces trois-là exactement.

- `reverse_2m` : pour 2m > 64 c'était une **boucle O(2m)** — 94 itérations de shift/masque/or sur
  256 bits, une fois par super-k-mer produit. Remplacée par une réversion **par limbes** (chaque
  mot de 64 bits passé dans le réseau d'échange existant, ordre des limbes inversé, puis
  réalignement) : O(#limbes), 2 réversions au lieu de 94 itérations.
- `phi` : mixage en `__uint128_t` quand 2m ≤ 128.
- `mmer_repeats` : roulement du m-mer en `__uint128_t`.

**Identité — et un piège trouvé par le balayage.** Les trois sont identiques *en valeur* (les bits
survivants sont les 2m bits bas, et un produit `__uint128_t` donne exactement les 128 bits bas).
Mais la première version de `reverse_2m` a **échoué le balayage k/m à k=40, m=33** : là,
2m = 66 alors que `2(2k−m) = 94` sélectionne une largeur de génération de **64 bits** — une
fenêtre de minimizer plus large que l'entier de travail. La marche par limbes décalait alors un
`uint64_t` de 64 (UB, résultat différent de la boucle). Le garde `w <= bits(kuint)` conserve la
boucle historique pour ces configurations dégénérées.

| Mesure (chr21, N=7) | Gain |
|---|---|
| construct k=95 t=1 | **+60,86 %** (10,83 s → 6,73 s) |
| construct k=95 t=6 | **+56,04 %** (2,74 s → 1,75 s) |
| query k=95 t=1 | **+49,29 %** |

Attribution par profil avant/après : `reverse_2m` ≈ −1,6 s, `mmer_repeats` ≈ −1,5 s,
`phi` ≈ −1,0 s (total −4,1 s, cohérent avec le mesuré). Aucun effet aux largeurs plus étroites
(branches `if constexpr (sizeof(kuint) > 16)`).

### L7 — court-circuit XOR dès 64 bits (`e30d64b`)

`QUERY_XOR_MIN_STORE_BYTES` valait 16, ce qui laissait les enregistrements `uint64` (k ≈ 22…40 au
bucketing par défaut — la configuration la plus courante) sur la boucle d'origine, qui remasque
**les deux** opérandes à chaque sonde. Le commentaire du code enregistrait « −0.0 % at k31 » ; la
re-mesure de zéro donne l'inverse. Le seuil est descendu à 8 : à 4 octets le gain est négatif
(−0,58 %), un pair 32 bits étant déjà une comparaison unique et bon marché.

| Mesure (N=9) | Gain |
|---|---|
| query chr1 k=31, p=0 / p=50 / p=100 | +3,26 % / +2,56 % / +2,29 % |
| query chr1 k=31, jeu 120 k reads | +3,64 % |
| query chr21 k=31 | +3,77 % |
| query chr21 k=15 (`uint32`, chemin inchangé) | +0,39 % |
| query chr1 k=31 t=6 | +1,53 % (non concluant à t=6) |

### L18 — réutilisation des buffers d'overlaps et de chaîne (`a77e22c`)

Trois sources d'allocation dans un chemin appelé (k−m) fois par bucket, sur des milliers de
buckets : `candidate_overlaps` était un **local** de la fonction (donc reparti de zéro à chaque
bucket) ; `valid_overlaps` était **redéclaré à l'intérieur de la boucle de colonnes** (masquant
une déclaration externe jamais utilisée) ; et `greedy_chaining` allouait `parent`, `tail` et son
résultat à chaque appel. Les deux vecteurs deviennent des membres et `greedy_chaining` gagne une
forme `greedy_chaining_into(begin, end, out, GreedyScratch&)`.

| Mesure (N=5–7) | Gain |
|---|---|
| setop union chr21 k=31 | +3,34 % |
| setop combiné (3 sorties) k=31 | +3,55 % |
| setop union chr1 k=31 | +1,89 % |
| setop union chr21 k=63 | +0,56 % (dilué : moins d'enregistrements, chaînage plus lourd) |

Coût : pic RSS setop chr1 −t6 68,9 → 72,5 MB (+5,2 %, ≈ 0,6 MB par worker), borné et constant.

### L19 — table de préfixe de minimizer par bucket (`76f6f17`)

**Le levier le plus rentable de la campagne.** L'attribution par ligne mettait 16 % du wall de
requête sur **une seule ligne** — la lecture des tailles préfixe/suffixe d'un enregistrement
sondé. Ce n'est pas de l'arithmétique : c'est le **défaut de cache L3 sur l'enregistrement**. La
requête est memory-latency-bound et payait ~log₂(bucket) tels défauts par groupe de colonnes.

Observation qui débloque : **tout masque de colonne couvre l'intégralité du slot de minimizer**
(tous les k-mers d'un super-k-mer contiennent son minimizer), donc un enregistrement dont le
minimizer diffère de celui de la requête compare inégal **à toutes les colonnes**. Les payloads
étant triés par `m_pair`, c'est-à-dire minimizer-major, les enregistrements candidats forment
**une seule plage contiguë**. Une petite table par bucket sur les t bits de tête du slot de
minimizer stocké la localise en O(1), et les dichotomies par colonne démarrent dessus.

- construite en une passe O(n) au chargement du bucket, sous le même verrou, avant le
  `release`-store qui publie le bucket ;
- ~n/8 entrées `uint32` (t choisi pour que la table reste ~1/8 du bucket, plafonné à 12 bits et à
  la largeur du minimizer) ; buckets < 512 enregistrements : t=0, comportement d'origine ;
- **RAM seulement — aucun changement du format disque, aucune reconstruction d'index** ;
- une plage vide répond « absent partout » **sans une seule sonde** : d'où le gain maximal sur les
  requêtes absentes.

| Mesure (chr1, N=7) | p=0 | p=50 | p=100 |
|---|---|---|---|
| query k=31 | **+55,15 %** | +45,53 % | +38,61 % |
| query k=63 | +27,52 % | +25,52 % | +22,33 % |

Compteurs avant → après (chr1 k=31, t=1) : cycles 5,94 G → 3,75 G (−37 %) ; cache-misses
108 M → 41 M (−62 %) ; LLC-load-misses 10,6 M → 7,5 M (−29 %) ; branch-misses 55,4 M → 27,2 M
(−51 %). Coût : pic RSS 618 → 627 MB (+1,3 %) à k=31, 642 → 646 MB (+0,6 %) à k=63.

**Limite connue** : à k=15/m=7 le minimizer ne fait que 2m = 14 bits et b = 12 en consomme 12,
laissant `mini_bits = 2`. La table ne peut alors resserrer que d'un facteur 4 — d'où le gain nul
mesuré à k=15 (−0,82 %). C'est structurel : à petit m, le bucketing a déjà consommé le minimizer.

### L4 — parallélisation du parse de requête (`41ad75b`)

Une fois la recherche ~2,5× moins chère, **le producteur mono-thread est devenu le plafond** : il
pèse ~31 % du travail CPU à k=31, ce qui borne la requête à ~3× quel que soit `-t`. Mesuré sur la
baseline : `-t4` = 2,7× et **`-t6` plus lent que `-t4`**.

Le thread lecteur ne fait plus que découper l'entrée ; **chaque worker parse ET interroge** sa
propre tranche. Deux formes d'unité de travail, toutes deux reproduisant exactement le flux
séquentiel de super-k-mers :

- une suite d'**enregistrements entiers** (le cas courant : les reads sont bien plus courts que la
  cible de chunk) — énumérés intégralement comme le fait `FileSkmerator`, donc aucune couture ;
- une **sous-plage [a,b)** d'un enregistrement plus long que la cible : le worker énumère
  [a−marge, b+marge] et garde les super-k-mers dont **l'indice de création** tombe dans [a,b) —
  même pavage et même marge (4(2k−m)) que la phase 1 de `ParallelConstruct` ; le dernier chunk
  réclame aussi le vidage de fin de séquence, dont l'indice de création dépasse L.

L'ordre est préservé par construction : l'itérateur produit par indice de création croissant et
porte un seul indice de création à travers les morceaux d'un super-k-mer ambigu. **Contrairement à
`construct` — où la phase 2 trie, donc seul l'ENSEMBLE devait être préservé — la requête exige la
SÉQUENCE exacte**, ce que la porte ci-dessous vérifie directement.

| chr1 k=31, 120 k reads | baseline | HEAD | |
|---|---|---|---|
| −t1 | 6,60 s | 3,98 s | +66,0 % |
| −t2 | 5,76 s | 2,11 s | **+173,3 %** |
| −t4 | 2,23 s | 1,34 s | +67,1 % |
| −t6 | 1,59 s | 1,10 s | +43,8 % |

L'inversion `-t6 < -t4` de la baseline a disparu. (Le gain à −t1 vient de L19.)

### L8 — tri de colonne sur clés pré-extraites (`d3ec019`)

Le comparateur du tri par colonne appelait `kmer_compare`, qui masque **les deux** opérandes à
chaque comparaison et déréférence l'énumération **deux fois** — deux accès aléatoires dans un
tableau d'enregistrements de 24 à 72 octets, O(n log n) fois par colonne. Extraire `masked_kmer`
une fois par élément dans un tableau contigu (clé, id) transforme le tri en balayage.

**L'identité est garantie par construction, pas par chance** : l'élément i du tableau de clés
correspond à l'élément i de l'ancien tableau d'ids, et les deux comparateurs renvoient le même
booléen pour les éléments correspondants ; introsort exécute donc les mêmes comparaisons et les
mêmes échanges et produit la même permutation — **y compris la résolution des ex æquo**, ce qui
compte car `std::unique` garde le premier de chaque groupe et cet id devient le `last_id` d'un
virtual skmer.

| Mesure | Gain |
|---|---|
| construct chr21 k=15 t=1 | +4,91 % (le plus dense en enregistrements ⇒ le plus de tri) |
| construct chr21 k=63 t=1 | +2,80 % |
| construct chr21 k=31 t=1 | +1,36 % (sous le plancher) |
| construct chr1 k=31 t=6 | +2,79 % |

---

## 5. Bilan final : baseline `572dc15` → HEAD `d3ec019`

Protocole alterné, médiane, N=5–7, P-cores physiques.

### construct (chr21 sauf mention)

| k | t | baseline | HEAD | gain |
|---|---|---|---|---|
| 15 | 1 | 4,411 s | 4,128 s | +6,86 % |
| 31 | 1 | 3,295 s | 2,975 s | +10,76 % |
| 63 | 1 | 3,769 s | 3,455 s | +9,09 % |
| 95 | 1 | 9,710 s | 6,504 s | **+49,28 %** |
| 31 | 6 | 0,804 s | 0,745 s | +7,98 % |
| 95 | 6 | 2,345 s | 1,597 s | **+46,89 %** |
| 31 (chr1) | 6 | 5,154 s | 4,715 s | +9,29 % |

### query (stream, chr1 sauf mention)

| k | présence | t | baseline | HEAD | gain |
|---|---|---|---|---|---|
| 31 | 0 % | 1 | 2,554 s | 1,658 s | **+54,01 %** |
| 31 | 50 % | 1 | 2,395 s | 1,615 s | +48,33 % |
| 31 | 100 % | 1 | 2,258 s | 1,644 s | +37,30 % |
| 63 | 0 % | 1 | 1,873 s | 1,467 s | +27,63 % |
| 63 | 50 % | 1 | 1,762 s | 1,403 s | +25,53 % |
| 63 | 100 % | 1 | 1,667 s | 1,345 s | +23,99 % |
| 15 (chr21) | 50 % | 1 | 7,325 s | 7,386 s | −0,82 % |
| 95 (chr21) | 50 % | 1 | 1,391 s | 0,919 s | **+51,40 %** |
| 31, 120 k reads | mixte | 2 | 10,53 s | 3,865 s | **+172,45 %** |
| 31, 120 k reads | mixte | 4 | 2,875 s | 2,347 s | +22,52 % |
| 31, 120 k reads | mixte | 6 | 2,090 s | 1,766 s | +18,38 % |

### setop (chr21)

| opération | baseline | HEAD | gain |
|---|---|---|---|
| union k=15 J=0,5 | 2,582 s | 2,606 s | −0,91 % |
| union k=31 J=0,5 | 2,074 s | 2,014 s | +2,95 % |
| union k=63 J=0,5 | 2,221 s | 2,200 s | +0,95 % |
| union k=95 J=0,5 | 2,863 s | 2,804 s | +2,11 % |
| union k=31 J=0,1 | 2,994 s | 2,928 s | +2,25 % |
| union k=31 J=0,9 | 1,490 s | 1,440 s | +3,46 % |
| `--sizes` k=31 | 0,344 s | 0,342 s | +0,56 % |
| combiné 5 sorties k=31 | 4,727 s | 4,616 s | +2,42 % |
| union k=31 t=6 | 0,608 s | 0,586 s | +3,68 % |

**Les set-ops sont le workload le moins amélioré**, et c'est cohérent avec la carte du temps :
71,8 % de leur coût matérialisé est la recompaction, dont le cœur est un algorithme (jointure
d'overlaps + chaînage + merge de colonnes), pas du gaspillage identifiable. Les deux redondances
trouvées ont été traitées (L18 retenu, L16/L17 rejetés).

### Portes de correction (rejouées après chaque levier, et sur HEAD)

- **213/213** tests gtest, zéro échec ;
- `construct` **byte-identical** baseline == HEAD(−t1) == HEAD(−t8) sur **11 configurations (k,m)**
  couvrant les 4 largeurs et les cas dégénérés `2m > bits(kuint)` (k=40/m=33, k=48/m=33) ;
- `setop` mono-op (4 opérations) **et** combiné (5 fichiers de sortie + `--sizes`) byte-identical
  vs baseline et entre −t1/−t8, sur k=15/31/63/95 ;
- `query` byte-identical vs baseline pour **−t 1, 2, 4, 6, 8, 16** sur 4 largeurs × 3 taux de
  présence, **et** sur chr21 entier en requête (séquences de 39 Mb ⇒ chemin « gros chunk » avec
  marges) ;
- oracles KMC `tests/setop_verif.sh` et `tests/setop_multi_verif.sh` : *ALL CHECKS PASSED* ;
- oracle de requête direct : les **4 554 269** k-mers distincts d'ecoli (ensemble KMC) tous
  rapportés présents, **0 absent** ; 200 000 k-mers aléatoires ⇒ **0 présent**, sortie identique à
  la baseline.

---

## 6. Points demandant ton arbitrage

### 6.1 Bug préexistant : segfault de `construct` quand `2m − b ≥ 8·sizeof(gen)`

**Trouvé en balayant k/m, présent à l'identique dans la baseline — je ne l'ai pas corrigé.**

`make_prefix_bucketing` route par `mini >> shift` avec `shift = 2m − b`. Quand `shift` atteint ou
dépasse la largeur de l'entier de génération, le décalage est un UB ; sur x86 le compte est
masqué modulo la largeur, l'id de bucket sort des bornes et le processus segfaulte.

Condition exacte, vérifiée sur 21 configurations :

| k | m | 2m | 2(2k−m) | genW | 2m ≥ 8·genW + b ? | rc |
|---|---|---|---|---|---|---|
| 25 | 21 | 42 | 58 | 4 B | 42 < 44 → non | 0 |
| 25 | 22 | 44 | 56 | 4 B | 44 ≥ 44 → **oui** | **139** |
| 42 | 37 | 74 | 94 | 8 B | 74 < 76 → non | 0 |
| 42 | 38 | 76 | 92 | 8 B | 76 ≥ 76 → **oui** | **139** |
| 70 | 68 | 136 | 144 | 16 B | 136 < 140 → non | 0 |
| 95 | 93 | 186 | 194 | 16 B | 186 ≥ 140 → **oui** | **139** |

Repro : `sskm construct -k 42 -m 38 -f <fa> -o /tmp/x.sskm -t 1`.

Je ne l'ai pas corrigé parce que le bon correctif est un **choix sémantique qui t'appartient** :
(a) router à pleine largeur comme le fait déjà `route_minimizer` côté lecteur, ou (b) refuser la
configuration avec un message clair. Note connexe : dans ce régime `minimizer()` tronque déjà le
minimizer à 64 bits via `to_kuint()`, donc (a) demande de regarder la chaîne complète.

### 6.2 Packing des enregistrements (`VSKMER_6`) — tu l'avais autorisé, je ne l'ai pas livré

Tu m'avais donné le feu vert pour un bump de format. **Je ne l'ai pas pris, pour deux raisons
mesurées :**

1. **Techniquement bloqué en l'état.** `__attribute__((packed))` n'a aucun effet : clang émet
   *"not packing field 'm_pair' as it is non-POD for the purposes of layout"* (`-Wpacked-non-pod`)
   parce que `Skmer::pair` a des constructeurs utilisateur. Les tailles restent 16/24/48/72. Le
   livrer demanderait soit de rendre `pair` POD (elle a ~10 constructeurs, utilisés partout), soit
   d'introduire un **type d'enregistrement disque distinct** avec conversion à chaque
   chargement — ce qui ralentirait query et setop.
2. **Le bénéfice attendu a fondu.** L19 a divisé par ~2,5 le nombre d'enregistrements touchés par
   requête, et les set-ops sont bornés par la recompaction, pas par la taille des enregistrements.
   Le gain résiduel serait surtout la **taille d'index** : −25 % (u32), −16,7 % (u64), −25 %
   (u128), −5,6 % (kuint256) — au prix de casser la comparabilité de tous tes benchs antérieurs.

Si la métrique bits/k-mer compte plus que la vitesse, la variante vraiment payante n'est pas le
packing serré mais de **loger les tailles préfixe/suffixe dans les bits libres du pair** : à
k=31/m=15/b=12 le pair a 128 − 82 = 46 bits libres pour 14 bits de tailles, ce qui donnerait un
enregistrement de **16 octets au lieu de 24 (−33 %)**, aligné et de pas puissance de deux. Dis-moi
si tu veux que je l'instruise.

### 6.3 E/S temporaires de la phase 1 (5,4 % de construct, mesuré)

`strace -c` sur `construct chr21 k=31 -t1` : **0,175 s de syscalls sur 3,2 s** (18 925 `openat`,
14 803 `write`, 18 907 `close`). `SkmerBucketWriter::flush_bucket` ouvre, écrit et ferme le
fichier de bucket **à chaque vidage**, avec un tampon par bucket de 8 Ko (plancher). Les deux
sorties évidentes coûtent cher : un cache de fd est inefficace (ordre de vidage quasi aléatoire
sur 4096 buckets), et augmenter le budget d'écriture multiplie 32 Mo × n_threads de tampons.
La sortie propre serait de changer la disposition des temporaires (un fichier séquentiel par
worker + un index de segments par bucket, lu en `pread` scatter en phase 2). C'est un changement
structurel : je te le laisse en décision.

### 6.4 Requête à petit m

À k=15/m=7 aucun de mes leviers ne mord (−0,82 %) : le minimizer fait 14 bits, `--buckets 4096`
en consomme 12, il ne reste que 2 bits pour la table de préfixe. Si ce régime t'importe,
la piste est de **découpler le nombre de buckets de la quotientation** à petit m, ou d'indexer
la table sur les premiers nucléotides de flanc plutôt que sur le minimizer.

---

## 7. Ce que je n'ai pas couvert

- **chm13 (3,1 Gb)** : hors périmètre convenu. Toutes les confirmations « grande taille » sont sur
  chr1 (224 Mb, index de 620–650 MB, donc bien au-delà des 24 MB de L3).
- **Mesures froides systématiques** : une seule série froide (L15). Le reste est chaud, et c'est
  indiqué. Sur un index qui ne tient pas en page cache, le classement des leviers pourrait
  changer (L19 réduit les sondes, donc devrait s'améliorer, mais je ne l'ai pas mesuré).
- **t ≥ 8** : la machine n'a que 6 P-cores ; au-delà, on mesure des E-cores à 3,8 GHz et le bruit
  explose. Les gains multithread rapportés s'arrêtent à t=6.
- **Le balayage Jaccard complet** (7 points) : j'ai mesuré 3 points (0,1 / 0,5 / 0,9).
- **`--max-ram` (bucketing adaptatif)** : jamais exercé par cette campagne.
- Le plancher de bruit à t=6 (5 %) m'a empêché de conclure sur plusieurs gains multithread
  potentiellement réels mais inférieurs à ce seuil.

## 8. Reproduire

```bash
cmake -S . -B build-optim -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=clang++-18
cmake --build build-optim -j"$(nproc)"
# banc de mesure (protocole alterné) et portes de correction :
source benchmark/results/latest/optim/bench.sh     # measure / measure_ab / evict
bash   benchmark/results/latest/optim/gates.sh     # 5 familles de portes, ~10 min
```

`bench.sh` et `gates.sh` vivent sous `results/latest/` (git-ignoré) : ce sont des outils de
campagne, pas du harnais permanent. À promouvoir dans `benchmark/scripts/` s'ils te sont utiles.
