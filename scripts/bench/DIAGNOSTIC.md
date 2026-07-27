# Audit d'optimisation MeanFieldHom.jl / TensND.jl — diagnostic mesuré

État au 2026-07-27. Machine : znver3, Julia 1.12.6, `JULIA_NUM_THREADS=1`,
`opt=2`. Baseline enregistrée sur un worktree `git` du commit `0cf9fd5`
**sans instrumentation** ; plancher de bruit 1,5–2,5 % selon la campagne.

Tous les chiffres ci-dessous sont **mesurés**, jamais estimés. Quand une
mesure s'est révélée fausse, elle est signalée comme telle plutôt que
supprimée.

---

## 0. Point de départ et ce qu'il a révélé

`_Qnn_direct` (`src/Core/green_kernel.jl:98-110`) contracte un tenseur d'ordre 4
avec le noyau de Green par une boucle 6-profonde `for p,q,r,s2,α,β in 1:3`,
soit 729 itérations × ~12 flops, là où la contraction se factorise en
`U = (C·n̂)·ξ` puis `B = U·K⁻¹·Uᵀ` (~100 flops).

Deux faits, vérifiés par lecture et par `grep` exhaustif sur `src/`, `ext/`,
`test/`, `docs/` :

1. **La forme factorisée existe déjà deux fois dans le dépôt** —
   `Cracks/green_residue.jl:33-72` (`Tncon → V → M·Vᵀ`) et
   `Core/green_helpers.jl:96-119`. `_Qnn_direct` est la version naïve laissée
   derrière.
2. **`_Qnn_direct` n'est appelé nulle part**, ni `_acoustic_tensor`, son unique
   consommateur. C'est du code mort : à supprimer, pas à optimiser.

Le balayage systématique a trouvé la même classe de problème ailleurs — et
trois catégories plus lourdes que la boucle initiale.

---

## 1. Le poste dominant : `hill_tensor` / `cod_tensor` appelés deux fois

Chaîne d'appel vérifiée ligne à ligne :

- `Schemes/mori_tanaka.jl` `_mt_4` appelle, **par phase**,
  `_phase_dilute_concentration` *et* `_phase_stiffness_contribution` ;
- le premier fait `strain_strain_loc(geom, P_i, P₀_proj)` ;
- le second passe par `Core.stiffness_contribution` → `contribution.jl:33-41`
  → **`strain_strain_loc(incl, C₁, C₀)`, arguments identiques** ;
- `strain_strain_loc` (`localization.jl:72-82`) fait
  `P = hill_tensor(incl, C₀; kw...)` — l'objet cher.

Facteur 2 exact sur le poste dominant, à chaque évaluation MT et à chaque
itération SC. Trois autres instances du même motif :

| instance | fichier | coût dupliqué |
|---|---|---|
| fissures (MT) | `_phase_compliance_contribution` puis `Cracks/compliance.jl:168-175` | `cod_tensor` ×2 |
| SC / ASC | `self_consistent.jl:136-140` → `_phase_stress_strain_average` recalcule `A_raw` dès que `sym ≠ NoSymmetrize` et `P_i` non isotrope | `hill_tensor` × nb de bins |
| `LayeredSphere` | `strain_strain_loc` + `stiffness_contribution` + `_membrane_surface_stress` interne | récurrences Hervé-Zaoui ×3 |
| `LayeredSpheroid` | `conductivity_contribution` appelle déjà les deux localisations | récurrences confocales ×3 |

### Correctif et gains mesurés (palier 1, gate bit-à-bit)

Nouvelles génériques `Core.loc_and_stiffness` / `Core.loc_and_stress_average`,
dispatchées sur la **classe d'inclusion** (repli sûr sur `AbstractInclusion`,
chemin rapide sur `AbstractEllipsoidalInclusion`), plus des bundles dédiés pour
les fissures, `LayeredSphere` et `LayeredSpheroid`.

| cas | temps | allocations | travail (compteurs) |
|---|---|---|---|
| `schemes/mt.porous.oblate.isosym` | **−50,5 %** | −14,2 % | hill 2→1 |
| `schemes/mt.aniso_matrix` | **−49,9 %** | −49,9 % | hill 2→1, nœuds 210→105 |
| `schemes/mt.crack.penny.tri` | **−49,2 %** | −49,9 % | cod 2→1, nœuds 210→105 |
| `schemes/mt.crack.penny` | −26,3 % | −21,0 % | cod 2→1 |
| `schemes/mt.theta_binned_ti.n20` | **−18,6 %** | −7,2 % | **hill 40→20** |
| `schemes/mt.iso2.sphere` | −2,3 % | −5,4 % | hill 2→1 |

Le canal compteurs confirme que le gain vient du travail supprimé et non d'un
changement de comportement de la quadrature adaptative : le nombre de nœuds
suit exactement le nombre d'appels Hill (210→105).

**Gate** : 67/67 cas bit-à-bit identiques (sha256 du rendu `%.17g`), suite
complète verte (7142 pass, 0 échec).

### Contrepartie mesurée, à assumer

Sur les cas les moins chers, la déduplication **coûte** plus qu'elle ne
rapporte : le tuple de retour `(A, N)` est alloué et boxé (la chaîne de
localisation est inférée `Any` par construction, voir §5), alors que le solve
Hill économisé est analytique et quasi gratuit.

| cas | temps | allocations |
|---|---|---|
| `schemes/mt.conductivity.iso2` (Hill 2ᵉ ordre iso analytique) | **+11,3 %** (629 → 700 ns) | **+22,2 %** (432 → 528 o) |

Mesuré directement, 2000 échantillons, min et médiane concordantes. C'est le
seul cas de la suite où la balance est négative ; l'arbitrage
(−50 % sur les cas chers contre +70 ns sur le cas le moins cher) me paraît
clairement favorable, mais il est réel et documenté ici plutôt que passé sous
silence.

Les autres petites hausses (`asc.stiffness` +3,1 %, `sc.newton` +3,6 %,
`sc.porous.sphere.phi30` +5,3 %) portent sur des corps où le raccourci
`sym isa NoSymmetrize` ne faisait **déjà** qu'un seul solve : il n'y avait rien
à dédupliquer, seul le tuple s'ajoute.

---

## 2. Incohérence pré-existante du milieu de référence

Découverte en écrivant le palier 1, indépendante de la demande initiale.

| helper | ordre | référence réellement passée |
|---|---|---|
| `_phase_dilute_concentration` | 4 et 2 | `P₀_proj = _project_matrix(P₀, sym)` |
| `_phase_stiffness_contribution` | 4 | `P₀_proj` |
| `_phase_stiffness_contribution` | **2** | **`P₀` brut** |
| `_phase_compliance_contribution` | **4 et 2** | **`P₀` brut** |

Donc dans un même schéma MT, `A_dil` et `N` d'une **même phase** voient deux
milieux de référence différents dès que `symmetrize ≠ NoSymmetrize`.

Invisible aux tests actuels : toute phase portant un `symmetrize` dans `test/`
a une matrice **et** une propriété de phase isotropes, donc
`_project_matrix(P₀, IsoSymmetrize())` est un no-op à ~1e-16 près.

**État** : le palier 1 préserve le comportement actuel — les bundles concernés
se replient sur les deux appels séparés quand `P₀_proj !== P₀` (le garde est
l'identité d'objet, ce qui couvre `NoSymmetrize` **et**
`TISymmetrize(matrix_projection = :none)`). L'uniformisation est le palier 2,
non encore appliqué.

---

## 3. Redondances dans les intégrandes (non encore corrigées)

| site | constat | correctif |
|---|---|---|
| `Core/green_helpers.jl:71-80` | `Kns[i,j] = Σ C_{ikjl}(n̂_k ξ_l + ξ_k n̂_l)` en 81 itérations. Par symétrie majeure `C_{ikjl}=C_{jlik}`, le 1ᵉʳ terme **est** `Vs[i,j]` et le 2ᵈ `Vs[j,i]` : `Kns == Vs + Vsᵀ`, zéro flop. `Ks` symétrique (moitié des entrées calculées deux fois). | supprimer l'accumulateur `skn` ; `Ks` sur `i ≤ j` |
| `Elasticity/hill_3d_aniso_residue.jl:96-111` | la boucle `for α in 1:21` passe **le même `Q` et le même `z`** 21 fois ; chaque appel refait `derivative(Q)` jusqu'à `Q4` et recalcule `sqrt(1+z²)`, `log(z+√)` | cache par (φ, racine) |
| `hill_3d_aniso_decuhr.jl:56-66`, `hill_3d_cylinder_aniso.jl:44-56` | boucle 4-profonde pour un `K` **symétrique** puis `Matrix` + `inv` LU par nœud | `_sym3_inv_acoustic` existe déjà (`hill_3d_aniso_nestedquadgk.jl:47-68`) |
| `hill_2d_aniso.jl:67-74` | 16 appels `quadgk` séparés, chacun évaluant l'intégrande complète (16 composantes) pour n'en garder qu'une | un `quadgk` vectoriel — **16×** |
| `Cracks/green_residue.jl:66-92` | `Bpoly` rempli pour les 9 `(i,k)` puis symétrisé a posteriori ; `acc = acc + …` alloue un `Polynomial` neuf à chacun des 54 termes | boucle `i ≤ k`, buffer réutilisé |
| `Cracks/cod_numerical.jl:42-56` | `Tncon`/`A` (invariants en φ) recalculés à chaque nœud — les fonctions sœurs `_direct` les hissent déjà correctement | contexte pré-calculé |
| `Schemes/self_consistent.jl:337,341,367,379` | 3 évaluations gratuites du résidu SC complet (= un `hill_tensor` par phase) par itération de Newton | reporter `r_new` |

Repère de coût, mesuré : `hill_tensor` triaxial/triclinique en
`:nestedquadgk` alloue **103 Mo** et prend 83 ms pour **13 665 évaluations**
d'intégrande. En `:residues`, 6,8 ms et 105 nœuds.

---

## 4. TensND — `TensOrtho` est plus lent que le chemin générique

Le constat le plus net de la campagne, mesuré :

| cas | temps | allocations |
|---|---|---|
| `tensnd/dcontract.ortho_ortho` | **13 860 ns** | 9 136 o |
| `tensnd/dcontract.iso_ortho` | 10 900 ns | 8 992 o |
| `tensnd/dcontract.gen_gen` (tenseur **générique**) | **20 ns** | 928 o |
| `tensnd/dcontract.ti_ti` | 11 ns | 224 o |
| `tensnd/dcontract.iso_iso` | 3 ns | 80 o |

Contracter deux tenseurs **orthotropes structurés** est **~690× plus lent** que
contracter deux tenseurs génériques denses. Même hiérarchie sur l'accès :

| cas | temps | allocations |
|---|---|---|
| `tensnd/getindex.ortho` | **2 792 ns** | 2 304 o |
| `tensnd/getindex.ti` | 17,8 ns | 96 o |
| `tensnd/getindex.iso` | 2,9 ns | 48 o |
| `tensnd/collect.ortho` | **48 810 ns** | 60 448 o |
| `tensnd/collect.iso` | 491 ns | 784 o |

Causes identifiées :

1. `Base.getindex(t::TensOrtho, i,j,k,l) = get_array(t)[i,j,k,l]`
   (`tens_walpole.jl:1025`) — une allocation `Array{T,4}` et une boucle de 81
   itérations à ~40 multiplications **par accès scalaire**. `TensOrtho` est le
   **seul type structuré resté sur le chemin dense** : `TensISO`,
   `TensTI{4,N=5/6/8}` et `TensTI{2,N=2/3}` ont tous reçu leur `getindex`
   fermé. Comme `TensOrtho <: AbstractArray`, tout parcours générique devient
   O(81²) pour une opération O(81).
2. **Aucun `dcontract` fermé pour `TensOrtho`** : tout `⊡` passe par
   `same_basis` → `change_tens` → `get_array` dense sur **les deux** opérandes
   avant la moindre arithmétique. Or l'en-tête du fichier
   (`tens_walpole.jl:863-870`) écrit déjà la structure : dans le repère
   matériau la matrice KM est bloc-diagonale
   `[3×3 sym] ⊕ diag(2C₄₄,2C₅₅,2C₆₆)`, donc `A⊡B` = un produit 3×3 plus
   3 produits scalaires. Jumeau structurel exact de l'`inv` déjà implémenté
   (`inv.ortho` : 8,6 ns — la preuve que la forme fermée marche).

Autres points TensND relevés (non mesurés individuellement) : 14 usages
d'`OMEinsum` dont 4 (`otimesu`, `otimesl`, `sotimes`) ne font **aucune
sommation** ; champs abstraits `Tens.basis::Basis`, `TensRotated.basis`,
`TensOrthogonal.basis`, `CoorSystemNum.{χ,R,Γ}_func::Function` ;
`best_sym_tens` matérialise `Array(get_array(t))` trois fois et résout deux
fois le même problème aux valeurs propres 3×3.

---

## 5. Bogues pré-existants trouvés en route

Chacun **vérifié sur un worktree du commit de référence** avant d'être affirmé.

1. **`Core._quadgk` est du code mort** (`Core/quadrature.jl:17`) — le wrapper
   que « tous les sous-modules devraient utiliser » d'après son propre
   commentaire d'en-tête. Tous appellent `QuadGK.quadgk` directement.
   Porte la liste du code mort à **15 fonctions**.
2. **`SelfConsistent` + phase `TensTI` sous `TISymmetrize`** →
   `MethodError: no method matching _hill_3d_ti_coaxial(::Ellipsoid{Spherical}, ::TensTI{4,Float64,8})`.
   L'estimation courante devient un TI à **8** paramètres (résultat de la
   symétrisation exacte) et le noyau analytique TI-coaxial n'a de méthodes que
   pour 5 et 6 paramètres.
3. **ForwardDiff à travers une propriété de phase `TensTI`** → le constructeur
   interne `TensTI{order,T,N}(::NTuple{N,T}, ::Tuple{T,T,T})` exige le **même**
   `T` pour les paramètres et pour l'axe ; or l'axe reste `Float64` quand les
   paramètres deviennent `Dual`. Bogue AD de TensND.
4. **Instabilité de type inhérente** : `hill_tensor` n'est pas inférable, par
   conception — `_resolve_algo(Val(method), incl, C₀)` résout l'algorithme à
   l'exécution, ce qui est précisément ce qui fait fonctionner `:auto` et le
   repli `NestedQuadGK` sous `Dual`. Toute la chaîne de localisation est donc
   `Any`. C'est la cause de la contrepartie du §1.

Les points 2 et 3 sont marqués `@test_broken` dans
`test/Schemes/test_loc_bundles.jl`, avec le diagnostic en commentaire.

---

## 6. Erreurs de méthode commises et corrigées

Consignées parce qu'elles conditionnent la confiance à accorder aux chiffres.

1. **Première instrumentation type-instable.** `_maybe_count(f)` renvoyait la
   closure de comptage depuis une branche : type de retour `Union`, donc
   dispatch dynamique **par nœud** — exactement le poste que les compteurs
   servent à mesurer. Remplacé par un `struct _CountingFn{F}` et un
   branchement une fois par appel `quadgk`. Vérifié ensuite : allocations
   **identiques à l'octet** (103 363 248 o sur `hill.nqgk.tri`).
2. **`@inferred` sur les bundles.** J'assertais une propriété que le code n'a
   jamais eue (§5.4). Remplacé par la propriété utile : le bundle alloue
   strictement moins que les deux appels séparés.
3. **`samples=7` dans le canal temps.** Sur des cas à forte pression GC de
   l'ordre de la dizaine de µs, une seule pause GC tombe dans les 7
   échantillons et gonfle le `minimum` d'un facteur 2. A produit deux
   fantômes : `asc.stiffness` **+103 %** (réel : +4,5 %) et
   `sc.porous.sphere.phi30` **−46 %** (réel : +3,6 %). Détectés parce que le
   canal compteurs montrait un travail **inchangé**. Canaux découplés :
   allocations à 7 échantillons (déterministes, `min == max` sur les 67 cas),
   temps jusqu'à 10 000 échantillons dans un budget de 2 s.
4. **Cas sub-µs non fiables.** Même après (3), `control/dilute_dual.iso2`
   affichait **+63 %** avec des allocations identiques à l'octet ; mesure
   directe à 2000 échantillons : **−1,7 %**. Mécanisme : `evals` choisi par
   `tune!` diffère entre campagnes (175 vs 180), et comparer un
   min-de-1-appel à un min-de-N-appels-amortis n'est pas comparer la même
   statistique. Le diff signale désormais `~evalsN→M` et **refuse** de
   déclarer un tel cas déplacé.

Sans le canal compteurs, j'aurais rapporté un doublement du coût de l'ASC.
Sans la mesure directe de recoupement, j'aurais rapporté une régression de
63 % sur un cas de contrôle que je n'ai pas touché.

---

## 7. Ce qui reste

| palier | contenu | état |
|---|---|---|
| P2 | uniformisation de la référence `P₀` (§2), commit séparé, gate 1e-14 | non commencé |
| P3 | noyaux de quadrature (§3) | non commencé |
| P4 | StaticArrays sur les chemins chauds | non commencé |
| P5 | TensND v0.2.6 (§4 et §5.3) | non commencé |
| P6 | suppression des 15 fonctions mortes (§0, §5.1) | non commencé |

Rien n'est commité.
