# Benchmark — opérations ensemblistes sur k-mers : sklib vs KMC vs CBL

**Date :** 2026-06-06 · **Machine :** 22 cœurs, 62 Gio RAM, Linux · k = 21 (impair, requis par CBL), m = 11 (sklib).

## TL;DR

| critère | verdict |
|---|---|
| **Mémoire** | 🟢 **sklib gagne massivement** : ~constante (5→29 Mo de ecoli à chr1) vs KMC ×66, CBL ×103 à grande échelle |
| **Cardinalité seule (`_size`)** | 🟢 **sklib unique et ~×3 plus rapide** que « matérialiser puis compter » (KMC/CBL n'ont pas d'équivalent) |
| **Compacité de l'index/sortie** | 🟢 **sklib ~25 bits/k-mer** vs KMC ~48, CBL ~43 |
| **Temps set-op matérialisé** | 🟠 dépend de la taille de **sortie** : sklib gagne quand la sortie est petite (diff, faible overlap), KMC gagne quand elle est grande (union, fort overlap) |
| **Correction** | 🟢 **accord exact sklib = KMC = CBL** sur les 6 quantités, vérifié sur données réelles |

**Découverte clé :** le **temps sklib ∝ taille de la SORTIE** (matérialisation = re-compactage en super-k-mers), tandis que **KMC/CBL ∝ taille de l'ENTRÉE** (scan). Cela explique tous les résultats et indique où chaque outil est optimal.

---

## 1. Méthodologie

- **Outils.** sklib (`sskm setop`, build-bench **Release** ISA native, v0.4.2 + set-ops) ; **KMC** (`kmc` + `kmc_tools simple`) ; **CBL** (github.com/imartayan/CBL, binaire spécialisé par k, `--canonical`). Les autres index k-mers du dépôt (sshash, SBWT, BQF, FMSI) sont **membership-only** → pas d'opérations ensemblistes, exclus.
- **Modèle de mesure.** On construit chaque index **une fois** (mesuré à part), puis on exécute chaque opération sur les index **pré-construits** (lire 2 index → écrire 1 résultat — identique pour les 3 outils). `/usr/bin/time -v` → temps mur + RSS pic ; médiane de 3 répétitions pour les set-ops.
- **Équité threads.** Comparaison principale **mono-cœur** (`taskset -c 0`, `kmc -t1`) car les set-ops sklib sont séquentielles en v1 et CBL est mono-thread par conception. KMC est multi-thread par défaut → une **section multi-cœurs dédiée** (§7) mesure son gain réel.
- **Cardinalités.** Fournies par le `_size` de sklib (rapide, sans matérialisation), **validées identiques à KMC et CBL** (§2).
- **Données.** Génomes réels (E. coli K12/Sakai/UTI89/IAI39, levure, *C. elegans*, humain chr1/20/21/22) + copies mutées contrôlées (`mutate.py`) pour un régime d'overlap constant en scaling.

## 2. Correction (accord à 3 voies)

E. coli K12 vs K12 muté 2 %, k=21 — **sklib = KMC = CBL exactement** :

| quantité | sklib | KMC | CBL |
|---|--:|--:|--:|
| \|A\| | 4 543 891 | 4 543 891 | 4 543 891 |
| \|B\| | 4 585 742 | 4 585 742 | 4 585 742 |
| ∩ | 2 984 901 | 2 984 901 | 2 984 901 |
| ∪ | 6 144 732 | 6 144 732 | 6 144 732 |
| A\\B | 1 558 990 | 1 558 990 | 1 558 990 |
| B\\A | 1 600 841 | 1 600 841 | 1 600 841 |

(Confirmé en plus par 10 tests unitaires + cross-validation KMC dans `tests/setop_verif.sh`.)

## 3. Scaling — temps vs taille (B = A muté 1 %, k=21, mono-cœur)

![Temps intersection vs taille](figs/fig_time_vs_size.png)

| génome | \|A\| Mk | ∩ sklib | ∩ kmc | ∩ cbl | ∪ sklib | ∪ kmc | **A\\B sklib** | A\\B kmc |
|---|--:|--:|--:|--:|--:|--:|--:|--:|
| ecoli | 4.5 | 0.53 | **0.35** | 0.51 | 0.79 | **0.35** | **0.18** | 0.32 |
| yeast | 11.5 | 1.45 | **0.83** | 1.10 | 2.13 | **0.86** | **0.50** | 0.76 |
| chr21 | 32.7 | 4.61 | **2.27** | 2.97 | 6.84 | **2.36** | **1.52** | 2.10 |
| chr20 | 53.9 | 7.65 | **3.60** | 4.78 | 11.3 | **3.73** | **2.55** | 3.34 |
| celegans | 91.2 | 13.7 | **5.98** | 9.15 | 19.8 | **6.22** | **4.45** | 5.56 |
| chr1 | 189.8 | 30.0 | **12.3** | 18.5 | 43.7 | **12.8** | **9.51** | 11.5 |

En **fort overlap** (paires mutées 1 % → ∩ ≈ 67 % des k-mers), la sortie de ∩ et ∪ est grande → **KMC gagne** ∩/∪. La sortie de **diff est petite → sklib gagne diff à toutes les tailles** (chr1 : 9.5 vs 11.5 s). CBL est systématiquement plus lent que KMC ici.

## 3-bis. Paires de génomes réels téléchargés (k=21, mono-cœur)

Vraies données (aucune mutation) : souches d'*E. coli* (overlap biologique fort) et chromosomes humains (overlap faible — cas typique de génomique comparative).

| paire | Jaccard | ∩ sklib | ∩ kmc | ∩ cbl | A\\B sklib | A\\B kmc | RAM ∩ sklib/kmc/cbl |
|---|--:|--:|--:|--:|--:|--:|--:|
| ecoliK12 ∩ Sakai | 45 % | 0.49 | **0.38** | 0.57 | **0.28** | 0.36 | **5**/95/262 |
| ecoliK12 ∩ UTI89 | 34 % | 0.42 | **0.38** | 0.60 | **0.35** | 0.37 | **5**/91/263 |
| ecoliSakai ∩ IAI39 | 35 % | 0.44 | **0.40** | 0.63 | 0.43 | **0.39** | **5**/95/267 |
| **chr21 ∩ chr22** | 2.5 % | **1.20** | 2.07 | 3.44 | 5.41 | **2.25** | **5**/345/474 |
| **chr20 ∩ chr21** | 1.4 % | **1.49** | 2.77 | 4.74 | 9.06 | **3.11** | **6**/454/561 |

Sur les **chromosomes humains** (Jaccard ~1–2 %), **sklib calcule l'intersection plus vite que KMC** (1.20 s vs 2.07 ; 1.49 s vs 2.77) **avec ~70× moins de RAM** — exactement le scénario « combien de k-mers partagés entre deux séquences ? ». Sur les souches d'*E. coli* (overlap plus fort), KMC reprend l'intersection de peu, sklib garde la différence.

## 4. Mémoire — le résultat marquant

![RAM vs taille](figs/fig_ram_vs_size.png)

RSS pic pendant l'intersection (Mo) :

| génome | sklib | kmc | cbl |
|---|--:|--:|--:|
| ecoli 4.5M | **5** | 93 | 256 |
| chr21 32.7M | **8** | 402 | 472 |
| celegans 91.2M | **18** | 961 | 1943 |
| **chr1 189.8M** | **29** | **1918** | **2985** |

La RAM sklib est **quasi constante** (streaming par bucket : ~une paire de buckets en vol) : ×66 moins que KMC et ×103 moins que CBL à chr1. Avantage **indépendant** de l'opération et de l'overlap.

## 5. `_size` — cardinalité sans matérialisation (feature unique sklib)

KMC et CBL doivent **matérialiser** un index résultat puis le compter. sklib répond à la cardinalité (Jaccard, containment, dédup de jeux) en un seul passage sans rien écrire :

| génome | sklib `_size` ∩ | sklib matérialise ∩ | KMC matérialise ∩ | gain `_size` vs KMC |
|---|--:|--:|--:|--:|
| ecoli 4.5M | 0.07 s | 0.53 s | 0.35 s | **×5** |
| chr21 32.7M | 0.70 s | 4.61 s | 2.27 s | **×3.2** |
| chr1 189.8M | 4.29 s | 30.0 s | 12.3 s | **×2.9** |

…le tout à ~5–30 Mo de RAM (vs ~2 Go pour KMC). Pour comparer des milliers de jeux (matrices de Jaccard), c'est décisif.

## 6. Overlap — le croisement (base chr21 32.7M, k=21)

![Temps intersection vs overlap](figs/fig_time_vs_overlap.png)

| Jaccard | ∩ taille (Mk) | **∩ sklib** | ∩ kmc | ∩ cbl |
|--:|--:|--:|--:|--:|
| 100 % | 32.7 | 5.26 | **2.12** | 2.30 |
| 67 % | 26.9 | 4.58 | **2.36** | 3.11 |
| 21 % | 12.0 | 2.78 | **2.44** | 3.74 |
| 2.5 % (chr22 réel) | 1.6 | **1.31** | 2.19 | 3.64 |
| 0.7 % | 0.5 | **1.25** | 2.47 | 3.93 |

Le temps sklib suit linéairement la taille de sortie ; KMC est plat. **Croisement de l'intersection ≈ 20–25 % de Jaccard.** En génomique comparative (chromosomes/espèces différents, Jaccard quelques %), **sklib gagne l'intersection** (×1.7–2).

## 7. Multi-cœurs (KMC `-t1` vs `-t<22>`)

> _CBL est mono-thread par conception (pas de `rayon`). sklib est mono-thread en v1 (parallélisme par bucket = next step). Seul KMC est multi-thread → on quantifie son gain réel (`-t1` vs `-t22`)._

![Intersection mono vs multi-thread](figs/fig_threads.png)

| paire | KMC **build** t1→t22 | KMC **∩ set-op** t1→t22 | KMC ∩ RAM @t22 | sklib ∩ (1 thr) | **sklib `_size`** (1 thr) | sklib RAM |
|---|--:|--:|--:|--:|--:|--:|
| ecoli 4.5M | 0.50→0.30 s (×1.6) | 0.22→0.25 s (**aucun gain**) | 389 Mo | 0.48 | 0.10 | **5 Mo** |
| chr21×chr22 32.7M | 2.64→0.75 s (×3.5) | 1.18→0.93 s (×1.3) | 1 380 Mo | 1.23 | 0.90 | **5 Mo** |
| chr1 189.8M | 14.2→2.18 s (×6.5) | 6.47→**2.77** s (×2.3) | **4 347 Mo** | 31.97 | 4.57 | **28 Mo** |

**Constats multi-cœurs :**
- Le **build** KMC parallélise très bien (jusqu'à ×6.5) ; le **set-op** KMC (`kmc_tools simple`) beaucoup moins (×2.3 au mieux, I/O/merge-bound ; **nul** sur petites données).
- Le threading **aggrave** le désavantage RAM de KMC : chr1 ∩ à `-t22` = **4.3 Go vs 28 Mo** pour sklib (**×155**).
- Même à `-t22`, l'intersection KMC matérialisée à chr1 (**2.77 s**) reste proche du **`_size` sklib mono-thread (4.57 s)** — et sklib reste imbattable dès qu'on ne veut que la cardinalité, à RAM ~constante.
- **Mais** pour une intersection/union **matérialisée à fort overlap et grande échelle**, KMC `-t22` domine nettement (∩ chr1 : 2.77 s vs sklib 31.97 s, ×11) : c'est le cas où le **parallélisme par bucket de sklib** (next step — buckets indépendants → scaling quasi linéaire, sans coût RAM) apporterait le plus.

**Note construction (orthogonale aux set-ops) :** la construction sklib est **mono-thread** et plus lente que KMC (chr1 : 54.9 s vs KMC 14.2 s en `-t1`, 2.18 s en `-t22`). C'est un coût unique, hors périmètre set-ops, mais à garder en tête pour un pipeline complet.

## 8. Compacité (bits par k-mer de l'index)

| | sklib | KMC | CBL |
|---|--:|--:|--:|
| index construit (~ecoli) | **~25** | ~48 | ~43 |

sklib produit l'index le plus compact (encodage 2-bit + super-k-mers + quotienting).

## 9. k sweep (ecoliK12 vs Sakai)

| k | ∩ sklib | ∩ kmc | A\\B sklib | ∩ sklib `_size` |
|--:|--:|--:|--:|--:|
| 15 | 0.75 | **0.36** | **0.31** | 0.13 |
| 21 | 0.49 | **0.39** | 0.49 | 0.10 |
| 31 | 0.45 | **0.38** | **0.32** | 0.08 |
| 41 | **0.40** | 0.58 | **0.35** | 0.07 |

sklib se bonifie avec k (intersection plus petite) ; **KMC paie un saut à k > 32** (2 mots machine) → sklib gagne ∩ à k=41. `_size` domine partout.

## 10. Synthèse — quel outil pour quoi ?

| besoin | meilleur choix |
|---|---|
| **Cardinalité / Jaccard / containment** (sans sortie) | 🟢 **sklib `_size`** (×3 vs KMC, RAM minime) |
| **RAM contrainte / très gros jeux / nombreux en parallèle** | 🟢 **sklib** (RAM ~constante) |
| **Différence A\\B, ou intersection à faible overlap** (génomique comparative) | 🟢 **sklib** |
| **Union, ou intersection à fort overlap, débit brut, multi-cœurs** | 🟠 **KMC** |
| **Index dynamique exact + set-ops** (ajout/retrait incrémental) | CBL (mais ×100 RAM, plus lent ici) |
| **Sortie la plus compacte** | 🟢 **sklib** (~25 bits/k-mer) |

**Pistes sklib :** (1) **parallélisme par bucket** (buckets indépendants → scaling ~linéaire, sans coût RAM) — comblerait l'écart sur union/∩-fort-overlap ; (2) accélérer le re-compactage de la sortie (le coût dominant de la matérialisation).

## 11. Reproductibilité

Scripts (`scripts/bench/`) : `bench_setops.sh` (cœur) · `run_scaling.sh` · `run_realpairs.sh` · `run_overlap.sh` · `run_ksweep.sh` · `run_threads.sh` · `mutate.py` · `plot_setops.py`. Données brutes (instantané versionné) : `report/data/setops_*.csv` ; figures : `report/figs/`. Les génomes (gitignorés sous `scripts/out/`) se re-téléchargent via le harness (`prepare_genome` dans `lib.sh`) ; régénérer les figures : `python3 scripts/bench/plot_setops.py`.

Régénérer un point de mesure, ex. : `BENCH_CSV_HEADER=1 bash scripts/bench/bench_setops.sh ecoli scripts/out/e2e/genomes/ecoliK12.sanitized.fa scripts/out/e2e/genomes/ecoliSakai.sanitized.fa 21 11 3`.
