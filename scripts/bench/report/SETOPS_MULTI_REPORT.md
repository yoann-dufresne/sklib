# Opérations ensemblistes combinées — une seule passe (sskm setop)

Une opération `setop` combinée calcule, en **une seule passe** de fusion sur les deux listes,
n'importe quel sous-ensemble de {A∩B, A∪B, A\B, B\A} — matérialisé et/ou simplement compté — au lieu
d'une passe par opération. Les trois relations disjointes (`both`/`only_a`/`only_b`) sortent déjà d'un
seul `merge_columns` ; la lecture des buckets, la construction du CSR par colonne et la fusion sont
**partagées**, et les quatre cardinalités sont gratuites (3 compteurs).

Ce rapport compare, sur génomes réels (E. coli → chromosome humain 1), le mode **combiné** au mode
**séquentiel** (lancer les opérations séparément), et à **KMC** — dont `kmc_tools simple` sait aussi
enchaîner plusieurs opérations en une passe, ce qui en fait le concurrent direct. FMSI et CBL n'ont
pas de mode combiné (séquentiel uniquement, RAM lourde) et ne sont pas inclus ici.

## Méthode

- Binaire release `build-bench/bin/sskm`, KMC/`kmc_tools` du PATH.
- Index construits **une fois** par paire (mesurés à part), puis opérations sur les index pré-construits.
- **V1 = monothread** : `taskset -c 0`, `-t 1` (sklib) et `-t1` (kmc). k=31, m=21 (records `__uint128_t`).
- Médiane de 2–3 réplicats (1 pour chr1). `result_kmers` = |A∩B|+|A∪B|+|A\B|+|B\A| (somme, autoritative
  via `sskm --sizes`, validée == KMC par `tests/setop_multi_verif.sh`).
- Données brutes : `scripts/bench/report/data/setops_multi_v1.csv` (générées par
  `THREADS="1" scripts/bench/run_setops_multi.sh`).

## Résultats V1 (monothread)

### Combiné vs séquentiel — facteur d'accélération

| paire | Σ k-mers | **`--sizes`** comb/séq | **matérialiser** comb/séq | KMC matérialiser comb/séq |
|---|--:|--:|--:|--:|
| ecoliK12 ∩ Sakai   | 14,3 M | 0,077 / 0,311 = **4,04×** | 1,45 / 1,62 = 1,11× | 0,51 / 1,53 = 3,02× |
| ecoliK12 ∩ UTI89   | 15,2 M | 0,083 / 0,338 = **4,07×** | 1,58 / 1,76 = 1,11× | 0,52 / 1,54 = 2,98× |
| ecoliSakai ∩ IAI39 | 16,0 M | 0,088 / 0,355 = **4,03×** | 1,63 / 1,83 = 1,13× | 0,55 / 1,63 = 2,94× |
| chr21 ∩ chr22      | 133,2 M | 0,592 / 2,354 = **3,98×** | 12,6 / 14,1 = 1,12× | 3,65 / 10,8 = 2,96× |
| chr20 ∩ chr21      | 182,5 M | 0,739 / 2,935 = **3,97×** | 17,1 / 19,0 = 1,11× | 4,92 / 14,8 = 3,02× |
| **chr1 ∩ chr1.mut1%** | 526,6 M | 2,55 / 10,34 = **4,05×** | 52,8 / 58,8 = 1,11× | 19,5 / 63,7 = 3,27× |

(temps en secondes ; comb = 1 passe combinée, séq = somme des 4 passes séparées.)

### Lecture

- **Comptage (`--sizes`) : ≈ 4× partout.** Quatre `*_size` séparés refont quatre fois la même fusion ;
  le combiné la fait **une fois**. Sans matérialisation, la fusion *est* tout le coût → gain net ≈ N.
  C'est le résultat marquant : obtenir les 4 cardinalités à l'échelle chr1 passe de 10,3 s à **2,55 s**.
- **Matérialisation : ≈ 1,11×** (compaction par défaut). La ré-compaction en super-k-mers est par
  sortie (~90 % du temps) et **ne se partage pas** : produire {∩, A\B, B\A} re-compacte exactement les
  mêmes partitions que séparément. Le combiné n'économise que les **3 fusions redondantes** sur 4 — une
  fusion vaut ≈ le temps `--sizes` (chr1 : ~2,5 s), d'où ~6 s gagnés sur 59 s.
- **Matérialisation `--no-compact` : ≈ 1,3×** (ecoliK12∩Sakai : 0,72 vs 0,95 s). Sans ré-compaction, la
  fusion partagée pèse plus lourd dans le total, donc le gain combiné monte — sans atteindre 4× car
  l'écriture des enregistrements reste par sortie.
- **KMC matérialiser : ≈ 3×.** Le coût d'une op KMC est surtout la relecture des deux bases ; le mode
  `simple` multi-op la partage, d'où ~3× — confirmant que « combiné > séquentiel » vaut pour les deux
  outils, l'ampleur dépendant de la part partageable.

### Mémoire et vitesse absolue (matérialisation combinée, mono)

| paire | temps sklib / KMC | RSS sklib / KMC | RAM ×moins (sklib) |
|---|--:|--:|--:|
| ecoliK12 ∩ Sakai | 1,45 s / 0,51 s | **7** / 157 Mo | **22×** |
| chr21 ∩ chr22    | 12,6 s / 3,65 s | **23** / 763 Mo | **33×** |
| chr20 ∩ chr21    | 17,1 s / 4,92 s | **28** / 943 Mo | **34×** |
| **chr1 ∩ chr1.mut1%** | 52,8 s / 19,5 s | **79** / 2706 Mo | **34×** |

En valeur absolue, la **matérialisation** de KMC est ~3× plus rapide que sklib (sklib paye l'assemblage
en super-k-mers), mais sklib tient en **~30× moins de RAM**. À l'inverse, pour le **comptage**, sklib
`--sizes` (2,55 s à l'échelle chr1, 13 Mo) n'a pas d'équivalent KMC : KMC doit matérialiser puis compter.

## Correction

`tests/setop_multi_verif.sh` valide qu'**un seul** appel combiné produit, octet de contenu pour octet,
les mêmes ensembles que les 4 opérations KMC — cardinalités (depuis `--sizes` et les compteurs par
fichier) **et** contenu — sur E. coli K12/Sakai et chr21/chr22 (66,6 M de k-mers d'union tous présents).
Les tests C++ prouvent en plus que le combiné est **octet pour octet identique** aux 4 fichiers mono
(tous nb de buckets et largeurs, dont `__uint128_t`) et indépendant du nombre de threads.

> Note : ce travail a débusqué et corrigé un bug de déterminisme **préexistant** — la struct `Skmer`
> écrivait 8 octets de *padding* non initialisés à la largeur `__uint128_t` (k≳32), rendant les fichiers
> non reproductibles octet à octet (y compris entre nombres de threads de construction). Corrigé sans
> changer la taille disque (commit `:bug: fix(skmer)`).

## Conclusion V1

Le mode combiné apporte un gain franc là où le travail est **partageable** : **≈ 4× pour le comptage**
multi-relations, à toutes les échelles jusqu'à chr1. Pour la **matérialisation**, le gain est modéré
(≈ 1,1× compacté, ≈ 1,3× `--no-compact`) car l'assemblage des super-k-mers reste par sortie ; le combiné
y supprime néanmoins les fusions redondantes. sklib conserve son avantage mémoire (~30× moins que KMC).
La V2 (multithread) mesure le passage à l'échelle par bucket de ces mêmes opérations.
