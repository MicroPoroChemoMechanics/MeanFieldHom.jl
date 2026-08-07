# Audit d'optimisation MeanFieldHomogenization.jl / TensND.jl — diagnostic mesuré

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

## 3. Redondances dans les intégrandes (palier 3, corrigé)

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

### Gains mesurés (palier 3, gate 1e-14)

| cas | temps | allocations | ce qui a changé |
|---|---|---|---|
| `kernels/hill.decuhr.tri.321` | **−62,0 %** | **−80,7 %** | `_sym3_inv_acoustic` au lieu de `Matrix` + `inv` LU par nœud |
| `kernels/hill2.aniso` | **−35,2 %** | −18,5 % | un `quadgk` vectoriel au lieu de 16 scalaires |
| `kernels/hill.dual.nqgk.tri` | **−21,5 %** | −0,0 % | `Kns = Vs + Vsᵀ` et `Ks` sur `i ≤ j` |
| `schemes/sc.newton` | −2,1 % | **−18,4 %** | report de `r_new`, `Tref` dérivé du 1ᵉʳ résidu de boucle |

Le compteur `nodes` est **identique** avant/après sur chaque cas de
quadrature (13 665, 21 825, 105, 315 …) : la quadrature adaptative n'a pas
changé de comportement, donc la comparaison de temps porte bien sur le même
travail.

Trois items du plan ont été mesurés comme non rentables ou trop invasifs et
n'ont **pas** été appliqués : le cache `prepare_logI`/`prepare_logz` de la
boucle 21 du chemin `:residues` (le plus invasif, gain attendu non vérifié),
la boucle `i ≤ k` de `Cracks/green_residue.jl`, et la forme à contexte
pré-calculé des `_Qnn_star_*`. Le chemin `:residues` sort du palier inchangé
(`hill.residues.tri.321` −0,4 %, dans le bruit).

---

## 3 bis. StaticArrays sur les chemins chauds (palier 4, gate 1e-14)

`_qnn_pair_components!` — la boucle la plus interne de tout le module
`Cracks` — écrivait dans un tampon appelant et construisait ~10 `Matrix{T}`
3×3 sur le tas **par nœud α** (`Vp, Vm, Kp, Km` par diffusion, `iKp, iKm`
via `_inv3`, puis quatre temporaires pour les deux `(V·iK)·Vᵀ`). Elle est
devenue une fonction **pure** renvoyant une `SMatrix{3,3,T}` ; `_A_and_Tn`,
`_phi_cache` et `_inv3` renvoient également du statique.

| cas | temps | allocations |
|---|---|---|
| `kernels/cod.nqgk.ellipse03.tri` | **−85,3 %** (4,35 ms → 639 µs) | **−99,4 %** (18,48 Mo → 111 Ko) |

C'est le plus gros gain unitaire de la campagne. Le checksum bouge de
**6,9e-16** — de la réassociation flottante pure, à 1 ULP, attendue dès qu'on
remplace un produit matriciel générique par sa forme statique déroulée.

Deux fausses pistes écartées en route, toutes deux dans le sens de la
lisibilité :

- une première version de `_A_and_Tn` en arithmétique d'indices `mod1`/`fld1`
  sur des tuples : illisible, et `ntuple(f, 27)` **sans `Val` est
  type-instable** — j'aurais introduit exactement la régression que je
  cherchais à supprimer. Remplacée par `MArray` → `SArray` avec les boucles
  d'origine intactes ;
- une closure `iK = (i,j) -> iKt[…]` dans la boucle chaude DECUHR — le piège
  de boxing déjà rencontré sur ce dépôt. Remplacée par de l'indexation
  directe de tuple.

Hors périmètre, assumé : les `zeros(T,3,3,3,3)` de queue de noyau (chemins
froids), `Viscoelasticity/` (tout est dimensionné par le nombre de pas de
temps), `LayeredSpheroids/` (troncature Legendre dynamique), et les
`Matrix{Polynomial{ComplexF64}}` (3×3 mais éléments non-`isbits` : une
`SMatrix` n'y supprimerait aucune allocation). Le portage `SVector{21,T}`
des retours d'intégrande des back-ends Hill n'a pas été fait — il n'était
pas nécessaire pour atteindre le gain, et le chemin `Integrals`/DECUHR exige
une vérification séparée du caractère mutable du tampon.

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

### Correctifs et gains mesurés (palier 5, TensND v0.2.6, gate 1e-14)

Deux changements seulement, et le second n'était pas dans le diagnostic
initial :

**(a) `getindex(::TensOrtho)` en forme fermée.** `_ortho_entry` est extraite
de `get_array` pour qu'un accès scalaire évalue **une** composante au lieu
des 81.

| cas | temps | allocations |
|---|---|---|
| `tensnd/getindex.ortho` | **−99,8 %** (2 792 → 4 ns) | **−95,8 %** |
| `tensnd/collect.ortho` | **−94,7 %** (48,8 → 2,47 µs) | **−98,6 %** |
| `tensnd/get_array.ortho` | −44,9 % | +0,0 % |

**(b) `tensor_or_array` était type-instable — le vrai gisement.** `dim` vient
de `size(tab, 1)`, donc c'est une valeur **d'exécution** : écrire
`Tensor{order, dim}(tab)` construit un type non concret à la compilation et
la construction devient entièrement dynamique — **3 147 ns et 3 120 o** pour
un tableau de 81 éléments, contre 311 ns pour produire ce tableau. Or tout
tenseur structuré rejoint la route générique par cette fonction
(`change_tens` → `same_basis` → chaque opération binaire), donc ce coût était
payé **deux fois par `⊡`** entre opérandes structurés. Le `dim` est
maintenant canalisé par `Val`.

| cas | temps | allocations |
|---|---|---|
| `tensnd/dcontract.iso_ortho` | **−58,2 %** (10,9 → 4,38 µs) | −46,3 % |
| `tensnd/dcontract.ortho_ortho` | **−57,0 %** (13,9 → 4,77 µs) | −45,5 % |

C'est un gain plus large que le cas `TensOrtho` initialement visé : il porte
sur **toute** paire d'opérandes qui retombe sur le chemin dense.

**(c) La forme fermée du `dcontract`, et le goulot qu'elle a révélé.**
Reprise dans un second temps, après avoir établi la convention au lieu de la
supposer.

`inv_KM` n'avait aucun problème de repère : `KM(t)` est la Kelvin-Mandel
**canonique**, `inv_KM` la relit comme telle, et le va-et-vient est exact
pour les huit formes (ordres 2 et 4, symétrique ou non, dim 2 et 3). Mon
prototype faisait `inv_KM(Mₐ·M_b)` avec `Mₐ`, `M_b` les KM du **repère
matériau** : il interprétait donc des composantes matériau comme canoniques.
En repère canonique `Q = I` et l'erreur était invisible (1,8e-12) ; en repère
tourné elle valait 6039. Il manquait la congruence.

Vérifiée sur trois repères plutôt que postulée :

    KM(t) == Q · KM_material(t) · Qᵀ          (5,3e-15)

où `Q` est la représentation Kelvin-Mandel de `R ⊠ˢ R`, c'est-à-dire
`KM(rot6(θ,ϕ,ψ))` — que TensND possédait déjà. `Q` est **orthogonale**
(2,2e-16) : c'est exactement ce que Kelvin-Mandel apporte sur Voigt, où la
matrice analogue ne l'est pas.

Le résultat n'est pas un `TensOrtho` : le produit de deux blocs 3×3
symétriques n'est pas symétrique sauf s'ils commutent, donc `A ⊡ B` est
orthotrope **sans symétrie majeure** — 12 constantes pour 9 stockées, écart
mesuré à 87,6, pas du bruit. Même élargissement que `TensTI{4}` N=5 → N=6.
La méthode renvoie donc le `TensCanonical` que la route générique produisait,
à 6,2e-16.

Et en profilant la forme fermée, les 1115 ns restants **n'étaient pas dans
l'algèbre** :

| | avant | après |
|---|---|---|
| `frame(A) == frame(B)` | **1117 ns** | 162 ns |
| tout le reste cumulé | ~130 ns | ~130 ns |

`AbstractBasis <: AbstractMatrix`, et son `getindex` passait par
`vecbasis(ℬ, :cov)` — la surcharge `Symbol`, qui construit `Val(var)` à
partir d'une valeur d'exécution : dispatch dynamique à **chaque accès
scalaire** (67 ns). Le `==` générique d'`AbstractArray` en faisait 18. Ce
coût était payé par tout parcours générique d'une base et par chaque
`_check_same_reference` — la garde coûtait plus cher que l'algèbre gardée.

| cas | avant v0.2.6 | après (b) | après (c) |
|---|---|---|---|
| `dcontract.ortho_ortho` | 13 860 ns / 9 136 o | 4 770 ns | **164 ns / 304 o** |
| `dcontract.iso_ortho` | 10 900 ns | 4 380 ns | **164 ns** |
| `inv_KM` 6×6 | 228 ns | — | **37,5 ns** |

Soit **−98,8 %** sur ce qui était désigné comme l'opportunité restante la
plus rentable. Le facteur ~690× face au tenseur générique tombe à ~3,3×.

Deux bogues `Dual` du §5 ont été corrigés au passage (constructeur `TensTI`
à eltypes mixtes, et `_ti8_to_ti6`) ; les deux `@test_broken` correspondants
sont devenus de vrais `@test`. Les points 3 à 7 du palier 5 planifié
(`best_sym_tens`, dé-einsum-ification, paramétrage concret des bases) n'ont
pas été abordés.

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

## 7. Campagne finale — récapitulatif

`--label=P7-ortho --baseline=baseline.json --gate=1e-14 --repeat-suite=2`,
67 cas, machine au repos.

```
24 déplacés, 0 échec de gate, 1 « régression » de contrôle,
21 non fiables (evals différents)
plancher de bruit (contrôles, p90 de |Δt|/t) = 0,8 %
```

**Gains** (au-delà du seuil « déplacé » = max(3×bruit, 3 %)) :

| cas | temps | allocations |
|---|---|---|
| `tensnd/getindex.ortho` | −99,7 % | −95,8 % |
| `tensnd/dcontract.ortho_ortho` | **−99,3 %** | **−94,9 %** |
| `tensnd/dcontract.iso_ortho` | **−99,3 %** | **−95,6 %** |
| `tensnd/collect.ortho` | −94,4 % | −98,6 % |
| `tensnd/inv_KM.gen` | **−87,0 %** | +0,0 % |
| `kernels/cod.nqgk.ellipse03.tri` | −85,3 % | **−99,4 %** |
| `kernels/hill.decuhr.tri.321` | −62,0 % | −80,7 % |
| `schemes/mt.aniso_matrix` | −50,8 % | −50,0 % |
| `schemes/mt.porous.oblate.isosym` | −50,2 % | −14,2 % |
| `schemes/mt.crack.penny.tri` | −49,7 % | −50,0 % |
| `schemes/mt.crack.penny` | −49,0 % | −35,3 % |
| `tensnd/get_array.ortho` | −44,4 % | +0,0 % |
| `kernels/hill2.aniso` | −35,2 % | −18,5 % |
| `kernels/hill.dual.nqgk.tri` | −21,5 % | −0,0 % |
| `schemes/mt.theta_binned_ti.n20` | −17,3 % | −7,2 % |

**Correction** : 63 cas sur 67 restent **bit-à-bit identiques** (`0,0e+00`).
Les quatre qui bougent le font tous par réassociation flottante —
`cod.nqgk.ellipse03.tri` à 6,9e-16, `hill.decuhr.tri.321` à 1,7e-18 (passage
en statique), `dcontract.iso_ortho` à 4,0e-17 et `dcontract.ortho_ortho` à
1,9e-17 (forme fermée) — tous très loin sous la tolérance 1e-14.

**La « régression » de contrôle n'en est pas une.** `control/alv.voigt.n50`
sort à **−5,6 %**, c'est-à-dire *plus rapide* ; le harnais signale tout écart
de contrôle sans regarder le signe. Vérifié plutôt que supposé : reproductible
sur cinq processus frais (−4,8 à −6,8 %), allocations identiques à l'octet,
checksum bit-à-bit, compteurs de travail inchangés. Un A/B annulant le seul
correctif de `bases.jl` montre que **ce n'est pas lui** ; le candidat restant
est `inv_KM`, que l'ALV appelle en boucle pour convertir ses blocs de Mandel.
Je ne l'ai pas isolé formellement.

C'est la limite du jeu de contrôles : il a été choisi en supposant que les
paliers ne toucheraient pas aux primitives partagées. `inv_KM`,
`tensor_or_array` et la comparaison de bases sont globales, donc un contrôle
peut légitimement bouger — dans le bon sens ici.

**Ce qui monte.** Le seul poste déterministe est
`schemes/mt.conductivity.iso2`, +22,2 % d'allocation — la contrepartie du
tuple de bundle décrite au §1, sur le cas le moins cher de la suite (528 o
au total). Le reste du cluster à +4/+6 % (`hill.nqgk.tri.321` +4,4 %,
`hill.nqgk.oblate.tri` +4,0 %, `sc.porous.sphere.phi30` +5,0 %,
`asc.stiffness` +4,5 %, `alv/trapezoidal.n50` +5,9 %) est du bruit machine :
le contrôle `differential.iso2` bouge de +2,5 % et `alv/trapezoidal.n50` ne
traverse que du code supprimé. J'ai vérifié le mécanisme suspecté plutôt que
de le supposer — `_counted_quadgk` s'infère au **même type concret** que
`QuadGK.quadgk` appelé directement, donc l'instrumentation n'introduit pas
d'instabilité sur ces chemins.

### Vérification de bout en bout

| | |
|---|---|
| suite MeanFieldHomogenization | **7154 / 7154** |
| suite TensND | verte (AD 63/63, projections TI/ORTHO 89/89, NLopt 46/46) |
| `benchmark_pichler` | **24 / 24** |
| `benchmark_hill_derivative` | **17 / 17** |
| `benchmark_nlayers` | §1-4, contraintes locales à 5,8e-16 |
| `benchmark_porous` | 134 / 140 — **identique au commit pré-campagne**, chiffre pour chiffre |
| build Documenter | exit 0, aucune docstring orpheline |

Les 6 échecs de `benchmark_porous` (DifferentialScheme, φ ≥ 0,50, erreur
relative croissant de 2,6e-03 à 6,2e-02) ont été rejoués sur un worktree du
commit `0cf9fd5` : **strictement les mêmes valeurs**. Écart pré-existant vis
à vis d'echoes, sans rapport avec cette campagne.

### Ce qui reste ouvert

| sujet | pourquoi ce n'est pas fait |
|---|---|
| conteneur orthotrope à 12 paramètres | `A ⊡ B` de deux `TensOrtho` n'a pas la symétrie majeure (§4c). Le résultat retombe donc sur `TensCanonical`, comme avant. Un conteneur dédié garderait la structure sur une chaîne de contractions, mais remonterait dans MFH — exactement le scénario du `MethodError` sur `TensTI{4,T,8}` (§5.2) |
| isoler le primitif derrière le `−5,6 %` d'`alv.voigt.n50` | vérifié réel et bénéfique, `bases.jl` écarté par A/B ; l'isolement exact n'a pas d'enjeu de risque |
| cache `prepare_logI`/`prepare_logz` du chemin `:residues` | le plus invasif des items du palier 3 ; le chemin sort de la campagne inchangé |
| `SVector{21,T}` pour les retours d'intégrande Hill | non nécessaire au gain obtenu ; le chemin `Integrals`/DECUHR demande une vérification séparée du tampon mutable |
| `best_sym_tens`, dé-einsum-ification, bases à type concret | paliers TensND 4 à 7, non abordés |
| tags `Dual` imbriqués à travers `NewtonDefault` | débloqué par le correctif `_sc_newton_seed`, révèle un problème suivant, précédemment inatteignable |
