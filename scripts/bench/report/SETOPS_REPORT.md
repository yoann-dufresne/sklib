# Benchmark — opérations ensemblistes sur k-mers : sklib vs KMC vs CBL

**Date :** 2026-06-07 · **Machine :** Intel Core Ultra 7 165H (22 cœurs), 62 Gio RAM · k = 21 (impair, requis par CBL), m = 11 (sklib). sklib **Release** v0.5.0, set-ops **~1.6× plus rapides** (commit `5540b36`).

## TL;DR

| critère | verdict |
|---|---|
| **Mémoire** | 🟢 **sklib gagne massivement** : ~constante (5→29 Mo de ecoli à chr1) vs KMC ×66, CBL ×103 à grande échelle |
| **Cardinalité seule (`_size`)** | 🟢 **sklib unique et ~×3 plus rapide** que « matérialiser puis compter » (KMC/CBL n'ont pas d'équivalent) |
| **Compacité de l'index/sortie** | 🟢 **sklib ~25 bits/k-mer** vs KMC ~48, CBL ~43 |
| **Temps set-op matérialisé** | 🟢 **~1.6× plus rapide qu'avant** : sklib gagne ∩ sur **toutes** les paires réelles (chr21∩chr22 : 0.81 vs KMC 2.04 s, ×2.5) et à overlap < ~40–50 %, et **diff partout** ; KMC ne garde que l'union et ∩ à très fort overlap |
| **Correction** | 🟢 **accord exact sklib = KMC = CBL** sur les 6 quantités, vérifié sur données réelles |

**Découverte clé :** le **temps sklib ∝ taille de la SORTIE** (matérialisation = re-compactage en super-k-mers), tandis que **KMC/CBL ∝ taille de l'ENTRÉE** (scan). Cela explique tous les résultats et indique où chaque outil est optimal.

---

## 1. Méthodologie

- **Outils.** sklib (`sskm setop`, build-bench **Release** ISA native, v0.4.2 + set-ops) ; **KMC** (`kmc` + `kmc_tools simple`) ; **CBL** (github.com/imartayan/CBL, binaire spécialisé par k, `--canonical`). **FMSI** *supporte aussi* les opérations ensemblistes (framework **f-MS expérimental** : `fmsi inter/union/diff`) — mesuré au **§9 (k-sweep)**, mais **20–100× plus lent** que sklib/KMC et bien plus gourmand en RAM (ex. ∩ ecoli : 11–13 s / 480 Mo vs sklib 0,3 s / 5 Mo) ; cardinalité validée vs KMC (union exacte ; ∩/diff < 0,1 %, écart de re-comptage). sshash, SBWT et BQF sont **membership-only** (pas d'opérations ensemblistes).
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
| ecoli | 4.5 | 0.43 | **0.36** | 0.53 | 0.51 | **0.36** | **0.12** | 0.33 |
| yeast | 11.5 | 0.95 | **0.82** | 1.10 | 1.35 | **0.86** | **0.33** | 0.75 |
| chr21 | 32.7 | 2.90 | **2.53** | 3.02 | 4.20 | **2.35** | **0.97** | 2.08 |
| chr20 | 53.9 | 4.91 | **3.60** | 4.76 | 7.07 | **3.74** | **1.66** | 3.36 |
| celegans | 91.2 | 8.54 | **5.93** | 9.15 | 12.31 | **6.17** | **2.89** | 5.52 |
| chr1 | 189.8 | 18.69 | **12.25** | 18.34 | 27.35 | **12.69** | **6.22** | 11.34 |

En **fort overlap** (paires mutées 1 % → ∩ ≈ 67 % des k-mers), la sortie de ∩/∪ est grande → **KMC garde ∩/∪**, mais après l'optimisation set-op (~1.6×) **l'écart s'est resserré** (chr1 ∩ : sklib **18.7** vs KMC 12.3 s, contre 30.0 avant). La sortie de **diff est petite → sklib gagne diff à toutes les tailles** (chr1 : 6.2 vs 11.3 s). CBL est systématiquement plus lent que KMC ici.

## 3-bis. Paires de génomes réels téléchargés (k=21, mono-cœur)

Vraies données (aucune mutation) : souches d'*E. coli* (overlap biologique fort) et chromosomes humains (overlap faible — cas typique de génomique comparative).

| paire | Jaccard | ∩ sklib | ∩ kmc | ∩ cbl | A\\B sklib | A\\B kmc | RAM ∩ sklib/kmc/cbl |
|---|--:|--:|--:|--:|--:|--:|--:|
| ecoliK12 ∩ Sakai | 45 % | **0.33** | 0.39 | 0.58 | **0.18** | 0.37 | **5**/95/262 |
| ecoliK12 ∩ UTI89 | 34 % | **0.28** | 0.37 | 0.61 | **0.24** | 0.37 | **5**/91/263 |
| ecoliSakai ∩ IAI39 | 35 % | **0.31** | 0.39 | 0.63 | **0.30** | 0.39 | **5**/95/267 |
| **chr21 ∩ chr22** | 2.5 % | **0.81** | 2.04 | 3.39 | 3.43 | **2.21** | **7**/345/474 |
| **chr20 ∩ chr21** | 1.4 % | **1.02** | 2.75 | 4.69 | 5.72 | **3.07** | **7**/454/561 |

Sur les **chromosomes humains** (Jaccard ~1–2 %), **sklib calcule l'intersection ~2,5× plus vite que KMC** (0.81 s vs 2.04 ; 1.02 s vs 2.75) **avec ~50–70× moins de RAM** — exactement le scénario « combien de k-mers partagés entre deux séquences ? ». Après l'optimisation set-op, sklib **gagne aussi l'intersection sur les souches d'*E. coli*** (fort overlap : 0.33 vs 0.39 s). KMC ne reprend l'avantage que sur la **différence à grande sortie** (chr20∩chr21 : 5.7 vs 3.1 s).

## 4. Mémoire — le résultat marquant

![RAM vs taille](figs/fig_ram_vs_size.png)

RSS pic pendant l'intersection (Mo) :

| génome | sklib | kmc | cbl |
|---|--:|--:|--:|
| ecoli 4.5M | **6** | 93 | 256 |
| chr21 32.7M | **10** | 402 | 472 |
| celegans 91.2M | **22** | 961 | 1943 |
| **chr1 189.8M** | **38** | **1918** | **2985** |

La RAM sklib est **quasi constante** (streaming par bucket : ~une paire de buckets en vol) : ×50 moins que KMC et ×79 moins que CBL à chr1. Avantage **indépendant** de l'opération et de l'overlap.

## 5. `_size` — cardinalité sans matérialisation (feature unique sklib)

KMC et CBL doivent **matérialiser** un index résultat puis le compter. sklib répond à la cardinalité (Jaccard, containment, dédup de jeux) en un seul passage sans rien écrire :

| génome | sklib `_size` ∩ | sklib matérialise ∩ | KMC matérialise ∩ | gain `_size` vs KMC |
|---|--:|--:|--:|--:|
| ecoli 4.5M | 0.06 s | 0.43 s | 0.36 s | **×6** |
| chr21 32.7M | 0.46 s | 2.90 s | 2.53 s | **×5.5** |
| chr1 189.8M | 2.91 s | 18.69 s | 12.25 s | **×4.2** |

…le tout à ~5–30 Mo de RAM (vs ~2 Go pour KMC). Pour comparer des milliers de jeux (matrices de Jaccard), c'est décisif.

## 6. Overlap — le croisement (base chr21 32.7M, k=21)

![Temps intersection vs overlap](figs/fig_time_vs_overlap.png)

| Jaccard | ∩ taille (Mk) | **∩ sklib** | ∩ kmc | ∩ cbl |
|--:|--:|--:|--:|--:|
| 100 % | 32.7 | 3.39 | **2.25** | 2.45 |
| ~83 % | 29.6 | 3.07 | **2.13** | 2.79 |
| ~70 % | 26.9 | 2.90 | **2.23** | 2.94 |
| ~22 % | 12.0 | **1.79** | 2.33 | 3.46 |
| 2.5 % (chr22 réel) | 1.6 | **0.82** | 2.04 | 3.41 |
| ~1 % | 0.5 | **0.78** | 2.34 | 3.74 |

Le temps sklib suit linéairement la taille de sortie ; KMC est plat. **Croisement de l'intersection ≈ 40–50 % de Jaccard** (après l'optimisation ~1.6× ; ~20–25 % auparavant). En génomique comparative (chromosomes/espèces différents, Jaccard quelques %), **sklib gagne largement l'intersection** (×2–3).

## 7. Multi-cœurs (KMC `-t1` vs `-t<22>`)

> _CBL est mono-thread par conception (pas de `rayon`). sklib est mono-thread en v1 (parallélisme par bucket = next step). Seul KMC est multi-thread → on quantifie son gain réel (`-t1` vs `-t22`)._

![Intersection mono vs multi-thread](figs/fig_threads.png)

| paire | KMC **build** t1→t22 | KMC **∩ set-op** t1→t22 | KMC ∩ RAM @t22 | sklib ∩ (1 thr) | **sklib `_size`** (1 thr) | sklib RAM |
|---|--:|--:|--:|--:|--:|--:|
| ecoli 4.5M | 0.52→0.33 s (×1.6) | 0.24→0.28 s (**aucun gain**) | 377 Mo | 0.38 | 0.09 | **6 Mo** |
| chr21×chr22 32.7M | 2.64→0.73 s (×3.6) | 1.14→0.98 s (×1.2) | 1 378 Mo | 0.84 | 0.63 | **6 Mo** |
| chr1 189.8M | 14.3→2.19 s (×6.5) | 6.51→**2.72** s (×2.4) | **4 349 Mo** | 19.49 | 2.98 | **38 Mo** |

**Constats multi-cœurs :**
- Le **build** KMC parallélise très bien (jusqu'à ×6.5) ; le **set-op** KMC (`kmc_tools simple`) beaucoup moins (×2.3 au mieux, I/O/merge-bound ; **nul** sur petites données).
- Le threading **aggrave** le désavantage RAM de KMC : chr1 ∩ à `-t22` = **4.3 Go vs 28 Mo** pour sklib (**×155**).
- Même à `-t22`, l'intersection KMC matérialisée à chr1 (**2.72 s**) reste proche du **`_size` sklib mono-thread (2.98 s)** — et sklib reste imbattable dès qu'on ne veut que la cardinalité, à RAM ~constante.
- **Mais** pour une intersection/union **matérialisée à fort overlap et grande échelle**, KMC `-t22` domine (∩ chr1 : 2.72 s vs sklib 19.49 s, ×7.2 — contre ×11 avant l'optimisation set-op) : c'est le cas où le **parallélisme par bucket de sklib** (next step — buckets indépendants → scaling quasi linéaire, sans coût RAM) apporterait le plus.

**Note construction (orthogonale aux set-ops) :** la construction sklib est **mono-thread** et plus lente que KMC (chr1 : 51.4 s vs KMC 14.3 s en `-t1`, 2.19 s en `-t22`) — mais à **RAM bien moindre** (138 Mo vs ~1.8 Go). C'est un coût unique, hors périmètre set-ops, mais à garder en tête pour un pipeline complet.

## 8. Compacité (bits par k-mer de l'index)

| | sklib | KMC | CBL |
|---|--:|--:|--:|
| index construit (~ecoli) | **~25** | ~48 | ~43 |

sklib produit l'index le plus compact (encodage 2-bit + super-k-mers + quotienting).

## 9. k sweep (ecoliK12 vs Sakai)

| k | ∩ sklib | ∩ kmc | ∩ cbl | ∩ fmsi | ∩ sklib `_size` |
|--:|--:|--:|--:|--:|--:|
| 15 | 0.43 | **0.33** | 0.55 | 11.0 | 0.10 |
| 21 | **0.33** | 0.39 | 0.63 | 11.6 | 0.08 |
| 31 | **0.31** | 0.41 | 0.74 | 11.0 | 0.07 |
| 41 | **0.29** | 0.57 | 0.74 | 4.6 | 0.07 |
| 59 | **0.33** | 0.53 | 0.81 | 13.8 | 0.08 |
| 63 | **0.29** | 0.52 | — | 12.3 | 0.07 |

sklib se bonifie avec k (intersection plus petite, **∩ < KMC dès k=21** après l'optimisation) ; **KMC paie un saut à k > 32**. À **grand k (59/63)** — où **CBL ne peut plus** (k impair ≤ 59) — sklib reste sub-seconde (∩ 0.29–0.33 s) et **FMSI**, seul autre outil set-op à ces k, est **20–40× plus lent** (∩ 11–14 s, RAM ~0.3–0.5 Go). `_size` domine partout (≤ 0.1 s).

## 10. Synthèse — quel outil pour quoi ?

| besoin | meilleur choix |
|---|---|
| **Cardinalité / Jaccard / containment** (sans sortie) | 🟢 **sklib `_size`** (×3 vs KMC, RAM minime) |
| **RAM contrainte / très gros jeux / nombreux en parallèle** | 🟢 **sklib** (RAM ~constante) |
| **Différence A\\B, ou intersection jusqu'à ~40–50 % d'overlap** (génomique comparative) | 🟢 **sklib** |
| **Union, ou intersection à fort overlap, débit brut, multi-cœurs** | 🟠 **KMC** |
| **Index dynamique exact + set-ops** (ajout/retrait incrémental) | CBL (mais ×100 RAM, plus lent ici) |
| **Sortie la plus compacte** | 🟢 **sklib** (~25 bits/k-mer) |

**Pistes sklib :** (1) **parallélisme par bucket** (buckets indépendants → scaling ~linéaire, sans coût RAM) — comblerait l'écart sur union/∩-fort-overlap ; (2) accélérer le re-compactage de la sortie (le coût dominant de la matérialisation).

## 11. Reproductibilité

Scripts (`scripts/bench/`) : `bench_setops.sh` (cœur) · `run_scaling.sh` · `run_realpairs.sh` · `run_overlap.sh` · `run_ksweep.sh` · `run_threads.sh` · `mutate.py` · `plot_setops.py`. Données brutes (instantané versionné) : `report/data/setops_*.csv` ; figures : `report/figs/`. Les génomes (gitignorés sous `scripts/out/`) se re-téléchargent via le harness (`prepare_genome` dans `lib.sh`) ; régénérer les figures : `python3 scripts/bench/plot_setops.py`.

Régénérer un point de mesure, ex. : `BENCH_CSV_HEADER=1 bash scripts/bench/bench_setops.sh ecoli scripts/out/e2e/genomes/ecoliK12.sanitized.fa scripts/out/e2e/genomes/ecoliSakai.sanitized.fa 21 11 3`.
