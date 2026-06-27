# Revue critique complète de MOOSE
### Mock Observation Of Synchrotron Emission — audit en vue de la diffusion publique

*Revue réalisée le 10 juin 2026 — base de code analysée : intégralité de `src/`, `test/`, `config/`, `python/`, `setup.jl`, `Project.toml`, `README.md`.*

**Convention de lecture** : chaque problème identifié est marqué d'une étiquette de sévérité —
🔴 **CRITIQUE** (résultats scientifiques faux), 🟠 **IMPORTANT** (robustesse/performance majeure), 🟡 **MOYEN** (qualité de code, maintenabilité), ⚪ **FAIBLE** (cosmétique). La section 9 consolide la priorisation.

---

## 1. Compréhension du projet

### 1.1 Objectif scientifique

MOOSE produit des **observations synthétiques d'émission synchrotron polarisée** à partir de cubes de simulations MHD (champ magnétique B1/B2/B3, densité, température, optionnellement densité électronique et fraction d'ionisation). Il calcule :

- les cubes Stokes **I (via Tb), Q, U** en fonction de la fréquence, avec ou sans rotation Faraday interne ;
- le cube de **mesure de rotation** RM(x, y, z) et ses cartes intégrées ;
- la **synthèse RM** (Brentjens & de Bruyn 2005) → cubes FDF (|F|, Re F, Im F) ;
- des produits dérivés : angle de polarisation, fraction de polarisation, DM/EM, Tb, cartes intégrées de B, spectres de puissance, moments statistiques ;
- un mode **HEALPix** pour les cartes plein ciel.

L'émissivité synchrotron suit Padovani et al. (2021) (tables `e_perp`/`e_para` pré-calculées, interpolées en (B⊥, ν)), la densité électronique suit soit Wolfire et al. (2003), soit ne ∝ x·nH, soit un cube externe, soit une constante.

### 1.2 Architecture actuelle

```
MOOSE.jl (module unique, ~35 include() à plat)
│
├── Point d'entrée ───────────────────────────────────────────────
│   run_moose() / run_moose_interactive()      (SyntheticObservations/MOOSE.jl)
│   MOOSEFromConfig.MOOSE_from_config(json)    (MOOSE_from_config.jl, sous-module)
│   MOOSE_cli.jl                               (parseur ARGS manuel, ~30 drapeaux)
│   python/moose_frontend.py                   (wrapper subprocess)
│
├── Cœur de traitement ───────────────────────────────────────────
│   _run_moose_processing(cfg::RunConfig)
│     └─ pour chaque simulation × ligne de visée :
│          ProcessSynchrotron (3 surcharges selon ne_option)
│            └─ _process_synchrotron_common
│                 ├─ ReadSimulation (FITS + permutation LOS)
│                 ├─ ne_builder (Wolfire / ∝nH / cube / constante)
│                 ├─ RM = 0.81·cumsum(ne·B∥·δl)
│                 ├─ QUnu3D / QUnuNoFaraday3D  (Q, U par canal)
│                 ├─ Tnu3D → BrightnessTemperature (I)
│                 ├─ Filter (passe-bande Fourier optionnel)
│                 ├─ _add_noise! (optionnel)
│                 ├─ RMSynthesis → FDF (optionnel)
│                 └─ écritures FITS (WriteDataOnDisk, Header)
│
├── Physique ─────────────────────────────────────────────────────
│   PhysicalParameters/ : Constants, RM, ElectronDensity, Bperp,
│   IntrinsicAngle, PolarizationAngle/Fraction, BrightnessTemperature,
│   ConversionJyK, pressure, Borientation, DM-EM …
│   Synchrotron/ : QUnu.jl (interpolateurs d'émissivité, boucles 3D),
│   Tnu.jl, ProcessSynchrotron.jl
│   Faraday/ : RMSynthesis.jl, FaradayParameters.jl (non branché)
│
├── E/S ──────────────────────────────────────────────────────────
│   FileIO/ : FITSUtils, ReadSimulation, Header, WriteDataOnDisk,
│   HealpixIO ; Frequencies/FreqFile.jl
│
├── Analyse ──────────────────────────────────────────────────────
│   Statistics/ : PowerSpectrum, Moments, RMS, EffectiveWidth,
│   Statistics.jl (SummarizeStats, non branché)
│   Filtering/Filter.jl
│
└── Support ──────────────────────────────────────────────────────
    Utils/ : ask_user, print_progress, MooseError, validations
    SyntheticObservations/ : DictHeaderParameters, InstrumentalParameters
```

### 1.3 Flux de données

1. **Configuration** : trois voies convergent vers la struct `RunConfig` (26 champs) — prompts interactifs, JSON (schéma plat ou imbriqué, validé par `build_config`), drapeaux CLI (fusionnés dans un JSON temporaire).
2. **Lecture** : `read_FITS_file` charge chaque cube en entier (`read(FITS(f)[1])`), puis `permute_dims` réoriente selon la LOS (axe 3 = ligne de visée). Les composantes B sont réassignées (B1, B2 = plan du ciel ; B3 = LOS).
3. **Physique par pixel** : pour chaque canal ν (tableau construit en **MHz**), l'émissivité (ε⊥ − ε∥ pour P, ε⊥ + ε∥ pour I) est interpolée sur (B⊥, ν), l'angle χ = ψ_src + RM·λ² appliqué, puis sommation le long de l'axe 3 → cartes Q(ν), U(ν), I(ν) converties en température de brillance.
4. **Post-traitement** : filtre interférométrique passe-bande (masque 0/1 en Fourier), bruit gaussien, synthèse RM (F(φ) = K Σ P exp(−2iφ(λ²−λ₀²))).
5. **Sorties** : FITS 3D/2D avec en-têtes construits depuis `DictHeaderParameters`, journal récapitulatif, copie de la config.

### 1.4 Dépendances fonctionnelles — points de fragilité structurels

- **Tout passe par des chaînes de caractères** : `ne_option ∈ {"1","2","3","4"}`, drapeaux `"Y"/"N"`, LOS `"x"/"y"/"z"`. Le dispatch de `ProcessSynchrotron` se fait par `if/elseif` sur des strings, pas par le système de types de Julia.
- **Couplage implicite par convention de noms de fichiers** : `ReadSimulation` attend `B1_*.fits`, `densityHp.fits`, etc. ; aucun mécanisme de découverte ni de message d'erreur dédié.
- **Code mort branché à moitié** : `FaradayParameters.rmsynthesis_parameters`, `getRMSF`, `SummarizeStats`, `WolfireConstants` interactif, plusieurs prompts d'`InstrumentalParameters` — écrits mais jamais appelés par le pipeline.
- **Module monolithique** : un seul namespace `MOOSE` avec 35 includes ; aucune frontière interne, toute fonction voit toutes les autres.

---

## 2. Détection des bugs

### 2.1 🔴 BUG-1 — Incohérence du mapping B1/B2 pour la LOS « y » (erreur de chiralité → signe de U faux)

**Où** : `src/FileIO/ReadSimulation.jl`, les deux méthodes `ReadSimulation`.

La méthode à 6 arguments applique pour `LOS == "y"` le mapping **cyclique** `(B1, B2, B3) = (Bz, Bx, By)`, tandis que la méthode à 5 arguments — **celle utilisée par le pipeline** — applique `(B1, B2, B3) = (Bx, Bz, By)`, qui est **anti-cyclique**.

**Pourquoi c'est faux** : passer d'un repère (x, y, z) à un repère « plan du ciel + LOS » doit préserver l'orientation (déterminant +1). Le mapping anti-cyclique inverse la chiralité du repère : l'angle intrinsèque ψ_src = atan(B2, B1) + π/2 est calculé dans un repère miroir.

**Conséquences** :
- Le **signe de Stokes U est inversé** pour toutes les observations en LOS y (Q est également affecté dès que ψ_src ≠ 0 mod π/2) ;
- les angles de polarisation, les FDF complexes (Re/Im), les gradients d'angle et toute statistique de polarisation croisée entre LOS x/y/z sont **incohérents entre eux** ;
- les deux méthodes donnant des résultats différents, tout utilisateur appelant l'une ou l'autre obtient une physique différente **silencieusement**.

**Correction** : centraliser le mapping dans une unique fonction, convention cyclique pour les trois LOS, et ne plus jamais le dupliquer :

```julia
"""
    los_basis(Bx, By, Bz, los) -> (B1, B2, B3)

Repère direct (plan du ciel ⊕ LOS), permutations cycliques :
los = "z" → (Bx, By, Bz) ; "x" → (By, Bz, Bx) ; "y" → (Bz, Bx, By).
"""
function los_basis(Bx, By, Bz, los::AbstractString)
    los == "z" && return (Bx, By, Bz)
    los == "x" && return (By, Bz, Bx)
    los == "y" && return (Bz, Bx, By)
    throw(MooseError("LOS inconnue : $los (attendu x, y ou z)"))
end
```

Les deux méthodes `ReadSimulation` appellent `los_basis`, et un **test de régression** vérifie que pour un champ B uniforme connu, ψ_src est identique (au signe de permutation près) sur les trois LOS. *Attention* : ce correctif change les résultats existants en LOS y — c'est voulu et doit être documenté dans le CHANGELOG (les résultats actuels sont erronés).

### 2.2 🔴 BUG-2 — Le « SNR » utilisateur est utilisé directement comme écart-type du bruit

**Où** : `src/Synchrotron/ProcessSynchrotron.jl`, `_add_noise!`.

```julia
noise = rand(Normal(0, Noise_nu), size(A))   # Noise_nu = valeur saisie comme "SNR"
```

**Pourquoi c'est faux** : l'utilisateur saisit un *rapport signal sur bruit* (interface et README parlent de SNR), mais la valeur est employée comme **σ absolu en Kelvin**. Demander SNR = 100 ajoute un bruit de σ = 100 K — c'est l'inverse de l'intention : plus le SNR demandé est grand, plus le bruit injecté est fort.

**Conséquences** : toutes les sorties « bruitées » publiées avec ce code ont un niveau de bruit sans rapport avec le SNR annoncé ; les études de détectabilité ou de robustesse statistique fondées dessus sont invalides.

**Correction** (préserve la fonctionnalité, change la sémantique vers celle documentée) :

```julia
function _add_noise!(A::AbstractArray, snr::Real; rng=Random.default_rng())
    snr > 0 || throw(MooseError("SNR doit être > 0, reçu $snr"))
    # σ défini canal par canal : signal de référence = rms du canal
    for k in axes(A, 3)
        ch = @view A[:, :, k]
        σ = sqrt(mean(abs2, ch)) / snr
        ch .+= σ .* randn(rng, eltype(ch), size(ch))
    end
    return A
end
```

Le choix du « signal de référence » (rms par canal, max de P, rms de P polarisé…) doit être **explicite dans l'en-tête FITS** (`NOISEREF`, `SNR`) et dans le journal. Ajouter aussi un `rng` injectable pour la reproductibilité (graine dans la config).

### 2.3 🔴 BUG-3 — `FFTW.fftfreq(n, Δx)` : la fréquence d'échantillonnage est inversée

**Où** : `src/Filtering/Filter.jl`, `instrument_bandpass_L`.

`fftfreq(n, fs)` attend la **fréquence d'échantillonnage** `fs = 1/Δx`, pas le pas `Δx`. Le code passe `Δx` directement : les fréquences spatiales générées sont fausses d'un facteur `Δx²`.

**Conséquences** : le masque passe-bande sélectionne des échelles qui n'ont rien à voir avec `Lcut_small`/`Lcut_large` demandés. Comme le masque reste un anneau 0/1, le résultat « ressemble » à un filtrage et l'erreur passe inaperçue. Toute publication de cartes filtrées avec des échelles annoncées en pc est fausse.

**Aggravation — confusion d'unités** : `_apply_synchrotron_filter!` transmet Δx en **pc**, alors que le README documente les coupures **en pixels**, et `Lcut_small = 1.0` est codé en dur. Trois conventions coexistent (pixels, pc, hard-codé).

**Correction** :

```julia
function instrument_bandpass_mask(nx, ny, Δx, Lmin, Lmax)
    kx = FFTW.fftfreq(nx, 1/Δx)        # fréquences spatiales correctes
    ky = FFTW.fftfreq(ny, 1/Δx)
    kmin, kmax = 1/Lmax, 1/Lmin        # Lmin < Lmax, mêmes unités que Δx
    [ (kmin ≤ hypot(kxi, kyj) ≤ kmax) for kxi in kx, kyj in ky ]
end
```

Décider d'une unité unique (recommandation : **pixels**, Δx = 1, avec conversion pc→pixels faite en amont et tracée dans le journal), exposer `Lcut_small` dans la config, et ajouter un test : un sinus pur de longueur d'onde L doit survivre si Lmin ≤ L ≤ Lmax et être annihilé sinon.

### 2.4 🔴 BUG-4 — `reshape` de la table d'émissivité avec dimensions inversées

**Où** : `src/Synchrotron/QUnu.jl` (`EmissivityInterpolator`) et `src/Synchrotron/Tnu.jl` (`TemperatureInterpolator`).

```julia
eps = reshape(df.e_perp .- df.e_para, (length(nu), length(B)))
itp = Spline2D(B, nu, eps)
```

`Spline2D(x, y, z)` attend `z[i, j] = f(x[i], y[j])`, soit une matrice `(length(B), length(nu))`. Le `reshape` produit la **transposée logique** : il ne fonctionne « par accident » que parce que la table fournie est carrée (même nombre de valeurs de B et de ν) **et** triée dans l'ordre qui compense l'inversion. De plus, le sens du tri du CSV (B-rapide vs ν-rapide) n'est jamais vérifié.

**Conséquences** : avec toute table non carrée → erreur de dimensions (au mieux) ; avec une table carrée mais ordonnée différemment → **émissivités silencieusement permutées**, soit des cartes I/Q/U entièrement fausses sans aucun avertissement.

**Correction** :

```julia
function emissivity_grid(df)
    B  = sort(unique(df.B));  ν = sort(unique(df.nu))
    nrow(df) == length(B)*length(ν) ||
        throw(MooseError("Table d'émissivité incomplète : $(nrow(df)) lignes ≠ $(length(B))×$(length(ν))"))
    perm = sortperm(collect(zip(df.nu, df.B)))   # ν externe, B interne → colonne-major (B, ν)
    eps  = reshape((df.e_perp .- df.e_para)[perm], (length(B), length(ν)))
    return B, ν, eps
end
```

avec vérification d'exhaustivité du produit cartésien et test unitaire sur une table rectangulaire synthétique f(B, ν) = B + 1000ν.

### 2.5 🔴 BUG-5 — Course de threads sur `_HEADER_PARAMS_CACHE`

**Où** : `src/FileIO/WriteDataOnDisk.jl`.

`WriteQUnu3D` lance deux `Threads.@spawn` (écriture Q et U) qui appellent chacun `get!(_HEADER_PARAMS_CACHE, key) do ... end` sur un `Dict` global **non protégé**. `Dict` n'est pas thread-safe : deux `get!` concurrents peuvent corrompre la table de hachage (résultats fantômes, `KeyError`, voire segfault selon la version de Julia).

**Conséquences** : crashs intermittents irreproductibles, ou en-têtes FITS mélangés entre produits — typiquement le genre de bug qui apparaît uniquement chez les utilisateurs avec `JULIA_NUM_THREADS` élevé.

**Correction** : soit pré-remplir le cache avant le `@spawn` (le plus simple, le dictionnaire devient en lecture seule), soit le protéger :

```julia
const _HEADER_LOCK = ReentrantLock()
header_params(key) = lock(_HEADER_LOCK) do
    get!(_build_header_params, _HEADER_PARAMS_CACHE, key)
end
```

### 2.6 🔴 BUG-6 — En-têtes FITS : unités fausses (corrompt l'interopérabilité)

**Où** : `src/SyntheticObservations/DictHeaderParameters.jl`, `src/FileIO/Header.jl`.

- Les cubes Qnu/Unu/Tnu déclarent `CUNIT3 = "Hz"` alors que l'axe spectral écrit (`nuArray`) est en **MHz** → tout outil standard (CASA, astropy, DS9) lira des fréquences fausses d'un facteur 10⁶ ;
- `"Pnumax"` a `BUNIT = "rad.m^{-2}"` alors que c'est une amplitude polarisée en **K** ;
- `intBLOS`/`intBtotal` : `"muG cm^2"` et `intBperp` : `"muG/cm^2"` — l'intégrale ∫B dl est en **µG·cm** (ou µG·pc) ;
- `"ne"` : `CTYPE3 = "VEL"` avec `CUNIT3 = "pc"` (contradictoire) ;
- `BLENGTH` n'est pas un mot-clé FITS standard (préférer `HIERARCH MOOSE BLENGTH` ou un commentaire) ;
- `CDELT3 = specarray[2] - specarray[1]` plante (`BoundsError`) si l'axe n'a qu'un canal.

**Conséquences** : les fichiers produits sont **scientifiquement ambigus pour quiconque d'autre que l'auteur** — rédhibitoire pour une diffusion publique.

**Correction** : convertir une fois pour toutes l'axe spectral en **Hz** dès sa construction (voir BUG-12), corriger chaque BUNIT/CTYPE, et garder `CDELT3 = length(a) ≥ 2 ? a[2]-a[1] : 0.0`. Ajouter un test qui relit chaque produit avec FITSIO et vérifie les mots-clés.

### 2.7 🟠 BUG-7 — Arguments fantômes : `PixelLength_pc`, `PixelLength_cm`, `DistanceArray`

**Où** : les trois surcharges de `ProcessSynchrotron`.

Ces arguments sont acceptés puis **ignorés** : les valeurs sont recalculées en interne à partir de `size(B1, 3)`. Le pré-calcul `los_pixel_scale(cfg.BoxLength_pc.x, cfg.BoxLength_pix.x)` dans `_run_moose_processing` est donc du code mort, et — plus grave — `BoxLength_pix` saisi par l'utilisateur **n'a aucun effet** : si le cube FITS n'a pas la taille annoncée, le pas en pc est silencieusement différent de celui demandé.

**Conséquences** : RM, DM, EM et Tb dépendent linéairement de δl ; une discordance cube/config passe inaperçue et fausse l'amplitude de tous les produits.

**Correction** : supprimer les arguments inutilisés, calculer δl à un seul endroit, et **valider** : `size(cube) == (BoxLength_pix...)` sinon erreur explicite. Gérer aussi les boîtes non cubiques (BoxLength_pc.x/y/z distincts) qui sont actuellement écrasées par un seul scalaire.

### 2.8 🟠 BUG-8 — `Tnu3D` : construction d'un interpolateur **par pixel et par fréquence**

**Où** : `src/Synchrotron/Tnu.jl`, `emissivity_total_at_frequency!`.

Un objet `linear_interpolation` est reconstruit à chaque appel dans la boucle interne. Pour un cube 512³ × 50 canaux, c'est ~6,7 milliards de constructions d'objets → des heures de GC pour un calcul qui devrait prendre des minutes. C'est aussi incohérent avec `QUnu.jl` qui, lui, pré-calcule un cache par fréquence (`build_emissivity_frequency_cache`).

**Correction** : répliquer exactement la stratégie de QUnu — pour chaque ν, extraire la coupe 1D ε(B) et interpoler avec une simple recherche binaire + interpolation linéaire en place (voir §4.3 pour la version optimisée).

### 2.9 🟠 BUG-9 — `read_FITS_file` : aucune validation, handle non fermé

**Où** : `src/FileIO/FITSUtils.jl`.

```julia
read_FITS_file(file) = read(FITS(file)[1])
```

- Le handle CFITSIO n'est fermé que par le GC → épuisement des descripteurs de fichiers sur les longues boucles multi-simulations ;
- aucun contrôle : HDU primaire vide (cas fréquent : données dans l'extension 1), cube 2D/4D, NaN, type inattendu → erreurs cryptiques bien plus loin dans le pipeline ;
- pas de message d'erreur indiquant *quel* fichier a échoué.

**Correction** :

```julia
function read_FITS_file(file::AbstractString; ndims_expected::Union{Int,Nothing}=3)
    isfile(file) || throw(MooseError("Fichier FITS introuvable : $file"))
    data = FITS(file, "r") do f
        hdu = findfirst(h -> h isa ImageHDU && length(size(h)) > 0, f)
        hdu === nothing && throw(MooseError("$file : aucun HDU image non vide"))
        read(f[hdu])
    end
    ndims_expected !== nothing && ndims(data) != ndims_expected &&
        throw(MooseError("$file : $(ndims(data))D lu, $(ndims_expected)D attendu"))
    any(!isfinite, data) && @warn "Valeurs non finies dans $file" count=count(!isfinite, data)
    return data
end
```

### 2.10 🟠 BUG-10 — Risques NaN/Inf non gardés dans la physique

- `PolarizationFraction` : `P ./ I` → **Inf/NaN** dès qu'un pixel a I = 0 (bords filtrés, masques). Garde : `ifelse(iszero(I), zero(P), P/I)` ou seuillage I > ε ;
- `Borientation` : `acos(arg)` avec `arg` issu d'un quotient flottant → `DomainError` pour `|arg| = 1 + 1e-16`. Garde : `acos(clamp(arg, -1, 1))` ;
- `Wolfire_ne` : `sqrt(Geff)` et `(T/100)^0.25` → NaN si Geff < 0 ou T < 0 (cubes de simulation avec valeurs négatives résiduelles près des chocs). Valider/clipper en entrée avec avertissement comptabilisé ;
- `_QUnu!` : si `RM` contient des NaN (propagé d'un ne ou B corrompu), `cos/sin(arg)` propage silencieusement → cartes Q/U trouées. Un `validate_finite(cube, name)` après chaque étape majeure coûte une passe O(N) négligeable et transforme des heures de débogage en un message clair.

### 2.11 🟠 BUG-11 — `MOOSEFromConfig.build_config` : `get(emiss_cfg, "path")` sans défaut

`get(dict, key)` à 2 arguments n'existe pas pour donner `nothing` — cela lève `MethodError`/`KeyError` selon l'usage. Tout JSON sans le champ `path` dans la section émissivité plante avec une erreur non explicite. Correction : `get(emiss_cfg, "path", nothing)` puis validation avec message dédié.

### 2.12 🟠 BUG-12 — L'axe des fréquences vit en MHz dans tout le pipeline

`nuArray` est construit en MHz, puis chaque consommateur reconvertit localement (`nui*1e6` dans `_QUnu!`, `nu_MHz*1e6` dans `BrightnessTemperature`, conversion implicite dans RMSynthesis qui, elle, attend des **Hz**…). Chaque conversion locale est une occasion d'oubli — c'est précisément l'origine de BUG-6 (CUNIT3). **Correction structurelle** : tout le pipeline interne en **SI (Hz)** ; les MHz n'existent que dans l'interface utilisateur (affichage, prompts). Idem pour pc/cm : une fonction `pc_to_cm` unique, pas de littéraux `3.086e18` dispersés.

### 2.13 🟡 BUG-13 — Convention d'intégration du RM non documentée + décalage d'une demi-cellule

`RM = cumsum(0.81 .* ne .* BLOS .* δl, dims=3)` :
- le sens du `cumsum` (indice 1 → n sur l'axe 3) définit *qui* est l'observateur ; ce n'est écrit nulle part, et la cohérence avec le sens de la permutation `permute_dims` n'est vérifiée par aucun test ;
- la cellule émettrice k subit la rotation de **sa propre cellule incluse** (cumsum jusqu'à k inclus) au lieu du milieu de cellule — biais systématique de ½ δRM par cellule, visible aux faibles résolutions.

**Correction** : documenter la convention (observateur côté k = 1 ou k = n), et utiliser le point milieu : `RM_k = cumsum_{j<k}(δRM_j) + δRM_k/2`.

### 2.14 🟡 BUG-14 — Divers

- `C_m = 2.99792458e8 # speed of light in cm.s^-1` (`Constants.jl`) : la **valeur** est en m/s, le commentaire dit cm/s. La valeur utilisée est correcte partout où je l'ai tracée (λ² en m², RM en rad/m²), mais ce commentaire est une bombe à retardement pour tout futur contributeur. Renommer `C_LIGHT_M_S` et corriger le commentaire.
- `CreateFreqFile` (`Frequencies/FreqFile.jl`) : le docstring n'est **pas attaché** à la fonction (chaîne libre suivie d'une ligne vide avant `function`) → invisible dans `?CreateFreqFile` ; `joinpath(repertory, ...)` après avoir déjà concaténé `* "/"` (double séparateur) ; pas de validation `num_freq ≥ 2`, `start < end`.
- `mv("T_nu.fits", "Tnu.fits")` dans `_process_synchrotron_common` : contournement fragile (échoue si le fichier existe déjà sans `force=true`, et révèle que le nom est mal construit en amont).
- `ConversionJyK` : argument `frame` accepté et ignoré ; unités GHz/arcsec implicites non documentées.
- `print_logo()` : les `sleep()` de l'animation s'exécutent aussi en mode batch/config → ralentit les jobs de cluster et pollue les logs SLURM de codes ANSI. Conditionner à `isatty(stdout)`.
- `_run_moose_processing` : dispatch sur `ne_option::String` par `if/elseif` sans branche `else` exhaustive → une valeur "5" passe silencieusement sans rien faire selon le chemin.
- Dépendances suspectes dans `Project.toml` : `Images`, `ImageFiltering`, `KernelDensity`, `StringEncodings` ne sont importées nulle part dans `src/` → alourdissent l'installation (temps de précompilation) pour rien.
- `getRMSF` et `rmsynthesis_parameters` existent mais ne sont jamais appelés par le pipeline : l'utilisateur n'a **aucun moyen automatique** de connaître la résolution en φ ni l'échelle maximale détectable de sa propre observation (voir §3.2).
- `Header.jl` : `CDELT3` plante si un seul canal (déjà cité en BUG-6) ; aucune écriture de `DATE`, `ORIGIN`, version de MOOSE — la provenance des fichiers est intraçable.

---

## 3. Audit scientifique

### 3.1 Synthèse RM (`Faraday/RMSynthesis.jl`)

**Ce qui est correct** : la forme F(φ) = K Σ_j P(λ²_j) exp(−2iφ(λ²_j − λ₀²)) avec K = 1/n_λ et λ₀² = moyenne des λ²_j est conforme à Brentjens & de Bruyn (2005) pour des poids uniformes. L'implémentation par produit matriciel (`mul!` sur P reshapé (nx·ny, n_λ)) est saine numériquement.

**Réserves** :
1. **Poids uniformes imposés**. B&dB définissent K = 1/Σw_j avec des poids w_j arbitraires (typiquement 1/σ²_j). Dès que du bruit hétérogène ou un drapeau de canaux (RFI) existe, l'absence de poids est une limitation scientifique réelle. → Ajouter `weights::AbstractVector = ones(nλ)`.
2. **λ₀² = mean(λ²)** : c'est le choix B&dB pour des poids uniformes (λ₀² = Σwλ²/Σw), donc cohérent — mais il doit être écrit dans l'en-tête FITS du cube FDF (`LAMBDA0SQ`), sinon personne ne peut interpréter les angles Re/Im.
3. **La dérotation s'arrête à λ₀²** : les angles de F(φ) sont référencés à λ₀², pas à λ = 0. C'est standard mais doit être documenté (l'angle intrinsèque s'obtient en dérotant de φ·λ₀²).
4. **`getRMSF` jamais branché** : le pipeline produit des FDF sans jamais écrire la **RMSF** correspondante ni sa FWHM (≈ 3.8/Δλ² — la constante 3.8 est la convention B&dB, OK). Une FDF sans RMSF est ininterprétable (déconvolution RM-CLEAN impossible). → Écrire systématiquement `RMSF.fits` + `FWHM_RMSF`, `PHI_MAX`, `MAXSCALE` dans l'en-tête, via `rmsynthesis_parameters` qui existe déjà.
5. **Échantillonnage en φ non validé** : l'utilisateur choisit `PhiArray` librement ; si δφ_grille > FWHM/3, les pics sont sous-résolus. `rmsynthesis_parameters` permettrait d'avertir automatiquement. De même, aucun avertissement si max|φ| demandé dépasse φ_max ≈ √3/δλ² (au-delà, la dépolarisation intra-canal rend F(φ) non fiable).
6. **δλ² calculé sur ν_min seul** dans `rmsynthesis_parameters` (`abs(2C²/ν³_min)·δν`) : c'est le canal le plus défavorable, donc conservateur — acceptable, à documenter.

### 3.2 RM et rotation Faraday interne

- Formule RM = 0.81 ∫ ne B∥ dl (cm⁻³, µG, pc) : prefacteur **correct**.
- `cumsum` : voir BUG-13 (convention de sens + demi-cellule).
- Dans `_QUnu!` : `χ = ψ_src + RM·(c/ν)²` avec ν convertie en Hz — **correct** (λ² en m², RM en rad/m²), mais dépend de la chaîne de conversions MHz (BUG-12).
- **Hypothèse physique implicite non documentée** : chaque cellule est traitée comme un émetteur derrière un écran de Faraday égal au cumsum — c'est le traitement standard du « Faraday interne » discretisé, mais il néglige la dépolarisation différentielle *intra-cellule* (facteur sinc(RM_cell·λ²)). Aux basses fréquences et fortes RM par cellule, cela surestime la polarisation. À documenter, et idéalement proposer en option le facteur de Burn intra-cellule.

### 3.3 Stokes I/Q/U et émissivité synchrotron

- Décomposition Q = Σ ε_P cos 2χ, U = Σ ε_P sin 2χ avec ε_P = ε⊥ − ε∥ et I = Σ(ε⊥ + ε∥) : **conforme** aux définitions standard (Padovani+2021).
- ψ_src = atan(B2, B1) + π/2 : l'orthogonalité E ⊥ B⊥ est correcte. **Mais** la validité dépend du repère (B1, B2) — d'où la criticité de BUG-1. Documenter aussi la convention IAU (angle compté depuis le Nord vers l'Est) et vérifier la cohérence avec l'orientation des axes FITS écrits.
- Température de brillance Tb = c²I/(2k_B ν²) (Rayleigh-Jeans, CGS) : **correcte** ; là encore ν doit être en Hz (conversion locale ×1e6, BUG-12).
- **Interpolation d'émissivité** : `Spline2D` (spline cubique Dierckx) sur (B⊥, ν) en **échelle linéaire** alors que ε(B, ν) varie en lois de puissance sur plusieurs décades. Une spline cubique en linéaire peut **osciller et produire des émissivités négatives** entre les nœuds aux bords de la table. → Interpoler en **log-log** (`log ε` sur `log B`, `log ν`) : la loi de puissance devient quasi linéaire, l'erreur d'interpolation chute de plusieurs ordres de grandeur et la positivité est garantie. C'est le changement scientifique le plus rentable hors bugs (à valider par comparaison avec l'intégration QuadGK directe sur quelques points).
- `linear_interp_extrapolated` : l'**extrapolation hors de la table** (B⊥ au-delà du max tabulé) est silencieuse. Une extrapolation linéaire d'une loi de puissance est fausse rapidement. → Compter et journaliser les pixels hors-table ; proposer `clamp` aux bornes par défaut + avertissement.
- B⊥ = √(B1² + B2²) (`Bperp.jl`) : correct.

### 3.4 Densité électronique

- **Wolfire et al. (2003)** : `ne = 2.4e-3·(ζ/1e-16)^½·(T/100)^¼·(G_eff)^½/ω_PAH + n·X_C` — forme conforme à l'éq. (presque universellement citée) de Wolfire 2003 pour le CNM/WNM. Réserves : (a) le **domaine de validité** (milieu neutre, T ≲ 10⁴ K) n'est pas vérifié — appliquée à des cellules de gaz chaud ionisé la formule est physiquement absurde ; ajouter un garde T_max documenté ou un avertissement ; (b) constantes ζ, G_eff, ω_PAH, X_C codées en dur dans le corps de la fonction au lieu d'être des paramètres de config (la version interactive `WolfireConstants` existe mais n'est pas branchée — incohérence) ; (c) NaN si T < 0 (BUG-10).
- ne ∝ x·nH (fraction d'ionisation) : correct si le cube x est bien une fraction ∈ [0, 1] — non validé.
- DM = ∫ne dl, EM = ∫ne² dl : formes correctes ; vérifier les unités d'en-tête (pc·cm⁻³, pc·cm⁻⁶).

### 3.5 Intégrations LOS, projections, cartes

- Sommation `sum(cube, dims=3) .* δl` : intégration en rectangle, cohérente avec le cumsum du RM. Acceptable ; documenter (pas de trapèze).
- `permute_dims` ([2,3,1] x, [3,1,2] y) : cohérent en soi, mais c'est le **pendant géométrique de BUG-1** — la permutation des axes spatiaux et celle des composantes vectorielles doivent être la *même* permutation. Test croisé indispensable : un dipôle placé en (i, j, k) doit se retrouver au pixel attendu pour chaque LOS.
- Cartes intégrées de B (`intBLOS` etc.) : ce sont des moyennes pondérées triviales — leurs unités d'en-tête sont fausses (BUG-6).

### 3.6 Statistiques, spectres, gradients

- `power_spectrum_2d`/`radial_psd` : périodogramme simple sans apodisation. Pour des cartes **non périodiques** (cubes filtrés, sous-régions), la fuite spectrale biaise la pente. → Option de fenêtre (Hann) + soustraction de la moyenne avant FFT ; documenter la normalisation (Parseval) qui n'est actuellement pas spécifiée.
- L'ajustement de pente est fait sur tout l'intervalle k sans exclusion des échelles bruitées/grille — exposer `kmin/kmax` d'ajustement.
- Étiquettes de tracé en français (« Spectre de puissance 2D », « pente ») : pour une diffusion publique, passer en anglais.
- `Moments.jl`, `RMS.jl`, `EffectiveWidth.jl` : implémentations standard, RAS, mais **aucun test**.
- `SummarizeStats` : trois blocs DataFrame copiés-collés, jamais appelé — supprimer ou factoriser et brancher.

### 3.7 Bilan de l'audit

| Domaine | Verdict |
|---|---|
| Synthèse RM (noyau) | ✅ Correct, mais RMSF/paramètres jamais fournis à l'utilisateur |
| RM interne | ✅ Prefacteur OK ; ⚠ convention cumsum + ½ cellule à corriger/documenter |
| Q/U/I | ⚠ Physique OK **sauf** LOS y (BUG-1) ; interpolation à passer en log-log |
| Tb, DM, EM | ✅ Correct (unités d'en-tête à corriger) |
| ne (Wolfire) | ⚠ Formule OK, domaine de validité non gardé |
| Filtrage | 🔴 Échelles fausses (BUG-3) |
| Bruit | 🔴 Sémantique SNR inversée (BUG-2) |
| Spectres | ⚠ Pas d'apodisation, normalisation non documentée |

---

## 4. Optimisation Julia

**Objectif fixé : cubes 512³ confortables, 2048³ possibles.** Ordres de grandeur : un cube Float32 512³ = 0,5 Gio ; 2048³ = 32 Gio. Le pipeline actuel charge simultanément B1, B2, B3, n, T (+ ne, x) **en Float64** (lecture FITS) : pour 2048³ c'est ~5 × 64 Gio = 320 Gio → impossible sur une station de travail. La stratégie doit donc combiner : (a) Float32 de bout en bout, (b) ne jamais matérialiser de cube temporaire évitable, (c) traitement par tranches pour les très grands cubes.

### 4.1 Empreinte mémoire — mesures structurelles

1. **Float32 par défaut**. Les cubes MHD sont en Float32 ; `read(FITS...)` les promeut souvent via les opérations Float64 ensuite (littéraux `0.81`, `2.0`…). Paramétrer le pipeline sur `T = eltype(cube)` et écrire les littéraux avec `T(0.81)` / `oftype`. Gain : ×2 mémoire partout, et ×2 de bande passante mémoire (le calcul est memory-bound).
2. **Éliminer les temporaires de `permutedims`**. `permutedims` **copie** le cube (0,5–32 Gio). Pour la LOS, utiliser `PermutedDimsArray` (vue sans copie) lorsque l'accès reste majoritairement séquentiel sur l'axe 3, ou ne permuter **qu'une fois** et libérer l'original (`cube = permutedims(cube, p)` puis l'ancien est GC-able ; aujourd'hui certaines fonctions gardent les deux vivants).
3. **Fusionner les broadcasts**. `deltaRM = RM_PREFACTOR .* ne .* BLOS .* PixelLength_pc` puis `cumsum(...)` matérialise deltaRM **et** RM (2 cubes). Remplacer par un cumsum en place fusionné :

```julia
function rm_cube!(RM::Array{T,3}, ne, B3, δl_pc) where T
    pref = T(0.81) * T(δl_pc)
    @inbounds for j in axes(RM,2), i in axes(RM,1)
        acc = zero(T)
        for k in axes(RM,3)
            δ = pref * ne[i,j,k] * B3[i,j,k]
            RM[i,j,k] = acc + δ/2      # point milieu (corrige BUG-13)
            acc += δ
        end
    end
    return RM
end
```

   Un seul cube de sortie, zéro temporaire, et le parcours k-interne est colonne-major-hostile → inverser la boucle (voir 4.2) ou accepter le coût (cumsum sur dims=3 a le même problème).
4. **Réutiliser les tampons Q/U/T entre canaux**. `QUnu3D` alloue les cartes par canal ; pré-allouer `Q[nx,ny,nν]` une fois et écrire dedans (déjà partiellement le cas) ; surtout, **ne pas garder** simultanément QUnu *et* QUnuNoFaraday si une seule version est demandée (aujourd'hui les deux chemins peuvent coexister).
5. **Mode « tranches » pour 2048³** (évolution majeure, voir §10) : toutes les opérations du pipeline (ne, δRM→cumsum, ε, sommation LOS) sont **séparables en colonnes (i, j)**. On peut donc streamer le cube par blocs de colonnes lus directement depuis FITS (`CFITSIO` permet la lecture par sous-régions) : empreinte mémoire O(nx·ny·n_bloc) au lieu de O(N³). C'est la seule voie réaliste vers 2048³ sur 64 Gio de RAM.

### 4.2 Boucles chaudes — `_QUnu!`

Version actuelle (simplifiée) : pour chaque canal, boucle threadée sur `CartesianIndices` avec `arg = 2(ψ + RM·f)` puis `cos(arg)`, `sin(arg)` séparés.

```julia
function QU_channel!(Qν, Uν, εP, ψ, RM, λ²::T) where T
    @inbounds Threads.@threads for j in axes(RM, 2)
        for i in axes(RM, 1)
            qacc = zero(T); uacc = zero(T)
            @simd for k in axes(RM, 3)
                s, c = sincos(2 * (ψ[i,j,k] + RM[i,j,k] * λ²))
                ε = εP[i,j,k]
                qacc = muladd(ε, c, qacc)
                uacc = muladd(ε, s, uacc)
            end
            Qν[i,j] = qacc; Uν[i,j] = uacc
        end
    end
end
```

Gains : `sincos` (1 appel au lieu de 2, ~×1,5 sur cette ligne) ; accumulation **dans la boucle LOS** → Q/U 2D écrits une fois (au lieu de matérialiser le cube cos/sin 3D par canal) ; `muladd` ; threading sur j (colonnes contiguës) plutôt que sur `CartesianIndices` (meilleure localité, pas de compteur atomique de progression dans la boucle chaude — le compteur atomique actuel sérialise partiellement les threads : le sortir au niveau canal).

**Point clé supplémentaire** : ψ_src et ε_P ne dépendent pas de φ ni — pour ψ — de ν : précalculer `2ψ` une fois (cube T), et ε_P par canal via le cache existant. Actuellement ψ est recalculé/relu tel quel mais `2.0*(...)` reste dans la boucle — micro, mais gratuit.

### 4.3 `Tnu3D` — correction de BUG-8 (le plus gros gain du code)

```julia
# Pré-calcul par canal : coupe 1D ε_tot(B) triée + recherche binaire inlinée
struct Eps1D{T}; B::Vector{T}; ε::Vector{T}; end
@inline function (e::Eps1D{T})(b) where T
    b ≤ e.B[1]   && return e.ε[1]          # clamp (cf. §3.3)
    b ≥ e.B[end] && return e.ε[end]
    hi = searchsortedfirst(e.B, b); lo = hi - 1
    w = (b - e.B[lo]) / (e.B[hi] - e.B[lo])
    return muladd(w, e.ε[hi] - e.ε[lo], e.ε[lo])
end
```

Boucle 3D identique à `QU_channel!` (accumulation LOS en registre). Sur 512³×50 canaux : on passe de milliards d'allocations à **zéro allocation dans la boucle**. Gain attendu : ×50–×500 sur Tnu.

### 4.4 Stabilité de type

- `RunConfig` : champs `Any`/`String` pour des drapeaux ("Y"/"N") et options ("1".."4") → chaque accès dans le code chaud est dynamique. Passer à `Bool` et `Symbol`/enum (`@enum NeOption wolfire propto_nH cube constant`) ; conversion une seule fois au parsing. Bonus : le compilateur peut spécialiser `ProcessSynchrotron` par dispatch.
- Fonctions retournant des tuples hétérogènes selon des branches string (ex. ReadSimulation qui retourne plus ou moins de cubes selon ne_option) : découper en méthodes dispatchées sur un type `NeModel` (voir §5).
- Vérifier avec `@code_warntype`/JET.jl les trois boucles chaudes (`_QUnu!`, `Tnu3D`, `RMSynthesis`) — la matrice de phaseurs de RMSynthesis est saine ; QUnu dépend du nettoyage RunConfig.

### 4.5 E/S FITS

- **Lire en mémoire-mappé/par région** pour les gros cubes (CFITSIO `read(hdu, :, :, k1:k2)`) — déjà cité en 4.1.
- **Écrire en Float32** : les sorties sont aujourd'hui écrites au type courant (souvent Float64) → fichiers ×2 trop gros sans gain scientifique. Paramètre `output_eltype` (défaut Float32).
- Les deux `Threads.@spawn` de `WriteQUnu3D` n'apportent rien (CFITSIO sérialise de toute façon les écritures disque) et créent BUG-5 → écrire séquentiellement.
- `apply_to_array_xy` (filtre) : FFT par tranche avec allocation à chaque tranche → `plan_fft` + tampon complexe réutilisé, ou `fft(A, (1,2))` en une passe sur le cube si la mémoire le permet.

### 4.6 RMSynthesis

L'implémentation `mul!` est déjà la bonne idée (BLAS). Améliorations : (a) construire la matrice de phaseurs en **ComplexF32** ; (b) traiter par blocs de pixels (tuiles de 10⁵ lignes) pour borner la mémoire de P reshapé ; (c) si n_φ et n_λ deviennent grands, l'algorithme NUFFT (φ ↔ λ² non uniforme) est la voie standard — à garder pour la vision long terme.

### 4.7 Divers

- `print_progress` appelé depuis les boucles threadées via compteur atomique : sortir l'affichage de la boucle chaude (rafraîchir au plus 10×/s depuis le thread principal).
- Constantes recalculées par canal (`(C_m/(ν·1e6))^2`) : précalculer `λ²::Vector` une fois — déjà fait dans RMSynthesis, pas dans QUnu3D.
- `cumsum(..., dims=3)` : sur l'axe 3, accès strided ; la version fusionnée de §4.1 le résout en même temps.
- Bench de référence à committer (`benchmark/`) : cube 256³ synthétique, mesures BenchmarkTools par étape, pour objectiver chaque optimisation et détecter les régressions.

---

## 5. Refactoring

### 5.1 Problèmes structurels

1. **Module monolithique** : 35 `include` à plat, tout est visible partout, aucune API délimitée. Conséquence directe : les deux `ReadSimulation` divergentes (BUG-1) ont pu coexister sans que rien ne le signale.
2. **Configuration tri-céphale** : interactive, JSON, CLI — trois chemins qui construisent `RunConfig` avec des validations différentes (le JSON est validé par `build_config`, l'interactif beaucoup moins). Un bug de validation se corrige donc trois fois.
3. **Sémantique stringly-typed** : "Y"/"N", "1".."4", "x/y/z" traversent tout le code.
4. **Fonctions à 10+ arguments positionnels** (`ProcessSynchrotron`, `_process_synchrotron_common`) dont certains ignorés (BUG-7) : signature illisible, erreurs d'ordre d'arguments indétectables.
5. **Duplication** : interpolateurs QUnu/Tnu quasi identiques ; trois surcharges ProcessSynchrotron à 80 % communes ; blocs DataFrame triplés dans SummarizeStats.

### 5.2 Architecture cible (préservant 100 % des fonctionnalités)

```
MOOSE.jl
├── Core/          types : RunConfig (typée), NeModel (types), LOSAxis (enum),
│                  Constants (SI), MooseError, validations
├── IO/            FITSUtils, ReadSimulation (+ los_basis unique), Header,
│                  WriteDataOnDisk, HealpixIO, FreqFile
├── Physics/       ElectronDensity, RM, Bfields, Angles, BrightnessTemperature,
│                  DMEM, JyK
├── Emissivity/    table (chargement/validation), interpolateurs log-log,
│                  caches par fréquence (fusion QUnu/Tnu)
├── Synchrotron/   QU_channel!, T_channel!, orchestration par canal
├── Faraday/       RMSynthesis, RMSF, rmsynthesis_parameters (branché)
├── PostProcess/   Filter, Noise, Statistics, PowerSpectrum, Moments
├── Pipeline/      run_pipeline(cfg) : l'unique orchestrateur
└── UI/            interactif, CLI, JSON → tous trois produisent un RunConfig
                   validé par UNE fonction validate(cfg)
```

Sous-modules Julia réels (`module IO ... end`) avec exports explicites : les frontières deviennent vérifiables et la documentation s'organise naturellement.

### 5.3 Refactorings ciblés

- **`NeModel` par dispatch** au lieu de `ne_option::String` :

```julia
abstract type NeModel end
struct WolfireNe   <: NeModel; ζ::Float64; Geff::Float64; ωPAH::Float64; XC::Float64; end
struct IonFracNe   <: NeModel end                 # ne = x · nH
struct CubeNe      <: NeModel; path::String; end
struct ConstantNe  <: NeModel; value::Float64; end

electron_density(m::WolfireNe, n, T) = @. 2.4e-3*sqrt(m.ζ/1e-16)*(T/100)^0.25*sqrt(m.Geff)/m.ωPAH + n*m.XC
electron_density(m::IonFracNe, n, x) = n .* x
# etc.
```

  Les trois surcharges de `ProcessSynchrotron` fusionnent en une seule fonction paramétrée par `m::NeModel` — suppression de ~200 lignes dupliquées et le `if/elseif` sans `else` (BUG-14) disparaît.
- **`RunConfig` typée** : `Bool` pour les drapeaux, `@enum`/types pour les options, `NamedTuple{(:x,:y,:z)}` conservés pour les boîtes, et `kwdef` avec défauts = ceux de `default_config.json` (une seule source de vérité ; générer le JSON de référence depuis la struct).
- **Arguments nommés** pour toute fonction > 4 paramètres.
- **Nommage** : mélange actuel CamelCase (`CreateFreqFile`, `ReadSimulation`) / snake_case (`run_moose`) / hongrois (`nuArray_MHz`). Convention Julia : `create_freq_file`, `read_simulation`… Garder des alias dépréciés (`@deprecate CreateFreqFile create_freq_file`) pour ne casser aucun script utilisateur.
- **Supprimer ou brancher le code mort** : `FaradayParameters` → brancher (§3.1) ; prompts inutilisés d'`InstrumentalParameters`, `WolfireConstants`, `SummarizeStats` → brancher ou supprimer ; dépendances fantômes du `Project.toml` → supprimer.
- **`MOOSE_cli.jl`** : remplacer le parseur manuel de 30 drapeaux par `Comonicon.jl` ou `ArgParse.jl` (aide auto-générée, types validés, complétion shell).

---

## 6. Expérience utilisateur (doctorants / postdocs / chercheurs)

### 6.1 Diagnostic

Le public visé veut : lancer vite une première observation synthétique, comprendre ce qui a été calculé, itérer sur les paramètres, scripter sur cluster. Aujourd'hui :

- le mode interactif pose ~25 questions séquentielles **sans valeurs par défaut affichées ni retour en arrière** — une faute de frappe à la question 24 oblige à tout recommencer ;
- les erreurs surviennent **tard** (fichier FITS manquant découvert après les prompts, parfois après de longs calculs) ;
- aucune **estimation préalable** mémoire/temps : l'utilisateur découvre l'OOM au bout d'une heure ;
- les paramètres dérivés essentiels (FWHM de la RMSF, φ_max, échelle max) ne sont jamais montrés (§3.1) ;
- la reprise d'une config existante est possible (JSON sauvegardé) mais non proposée par le mode interactif.

### 6.2 Recommandations

1. **Validation précoce et complète ("fail fast")** : à la fin de la saisie/du parsing, et *avant tout calcul* : existence et dimensions de tous les FITS, cohérence BoxLength_pix vs cubes, table d'émissivité couvrant [ν_min, ν_max] et la plage de B⊥ du cube, droits d'écriture sur le répertoire de sortie, estimation mémoire (`n_cubes × nx·ny·nz × 4 octets` + canaux) comparée à `Sys.free_memory()` avec avertissement explicite.
2. **Récapitulatif avant exécution** : tableau des paramètres + dérivés (δl, λ² min/max, FWHM RMSF, φ_max, mémoire estimée, nombre de tâches) et confirmation `[Entrée pour lancer / e pour éditer / q pour quitter]`.
3. **Défauts intelligents** : chaque prompt affiche son défaut (`Fréquence min [MHz] (défaut : 100) :`) ; Entrée = défaut ; le mode interactif propose d'emblée « ↵ recharger la dernière config / n nouvelle config ».
4. **Assistant de config = générateur de JSON** : le wizard interactif doit avoir pour *seul* effet de produire un JSON, ensuite exécuté par le même chemin que `MOOSE_from_config` — unification des trois voies (cf. §5.1-2) et reproductibilité gratuite.
5. **Presets** : `moose --preset lofar-hba`, `--preset askap-possum`… (gammes de fréquences réalistes) — trivial à fournir, très apprécié des doctorants.
6. **Messages d'erreur actionnables** : chaque `MooseError` doit dire *quoi faire* (« B1_sim42.fits introuvable dans /data — vérifiez `simulation_dir` ou le motif de nommage `B1_<nom>.fits` »).
7. **Graine aléatoire** dans la config (bruit reproductible) + version de MOOSE et hash de config dans tous les en-têtes FITS et le journal.

---

## 7. Interface utilisateur (CLI moderne)

Rester en CLI est le bon choix pour le public cluster. Améliorations concrètes :

1. **Menus numérotés robustes** : les `ask_user` actuels re-bouclent sur entrée invalide mais sans contexte. Standardiser un composant unique :

```
Densité électronique
  [1] Wolfire et al. (2003)   [2] ne ∝ x·nH
  [3] Cube externe            [4] Constante
Choix [1] :
```

   (défaut entre crochets, validation immédiate, `?` affiche une aide d'une ligne par option).
2. **Barre de progression avec ETA** : `print_progress` n'affiche pas d'ETA. Recommandation : `ProgressMeter.jl` (ETA, taux, compatible threads via `next!` thread-safe, se désactive proprement si `!isatty`). Une barre par étape (lecture, RM, Q/U, Tnu, filtrage, synthèse RM, écriture), pas une barre globale opaque.
3. **Journalisation propre** : remplacer les `println` épars par le système `Logging` (déjà partiellement utilisé : `@info` dans CreateFreqFile) avec deux sorties — console (niveau Info, concise, couleurs si TTY) et fichier `moose_run.log` (niveau Debug, horodaté). Les codes ANSI et `sleep` du logo désactivés hors TTY (BUG-14).
4. **Résumé final** : `write_summary_log` existe — l'enrichir : durée par étape (`format_duration` existe), liste des fichiers produits avec tailles, paramètres dérivés (FWHM RMSF…), avertissements accumulés (pixels hors table d'émissivité, valeurs non finies corrigées), et l'afficher aussi à l'écran en tableau compact.
5. **Codes de sortie** : le CLI a déjà des codes via MooseError — les documenter (`--help`) et garantir 0 = succès, 1 = erreur config, 2 = erreur E/S, 3 = erreur calcul.
6. **`moose --dry-run`** : exécute toute la validation + récapitulatif (§6.2) sans calcul — indispensable avant de soumettre un job SLURM de 10 h.

---

## 8. Documentation

### 8.1 État des lieux

- Docstrings inégales : certaines complètes (`CreateFreqFile` — mais détachée de sa fonction, BUG-14), beaucoup absentes (RMSynthesis, QUnu3D, ProcessSynchrotron — le cœur scientifique !) ;
- aucune doc des **conventions** (sens du cumsum RM, repère LOS, angle IAU, unités internes MHz/pc) — précisément les points où des bugs ont été trouvés ;
- README riche mais qui documente des comportements faux (coupures du filtre « en pixels » alors que le code passe des pc) ;
- pas de site Documenter.jl, pas de tutoriel exécutable.

### 8.2 Modèle de docstring (à appliquer au cœur scientifique en priorité)

```julia
"""
    RMSynthesis(Q, U, ν_Hz, φ) -> (absF, reF, imF)

Synthèse de mesure de rotation (Brentjens & de Bruyn 2005, éq. 25–38).

Calcule F(φ) = (1/n_ν) Σⱼ P(λ²ⱼ) exp[−2iφ(λ²ⱼ − λ₀²)] avec P = Q + iU
et λ₀² = ⟨λ²⟩ (poids uniformes).

# Arguments
- `Q`, `U` : cubes (nx, ny, nν) en K.
- `ν_Hz`   : fréquences des canaux **en Hz**, croissantes.
- `φ`      : grille de profondeur Faraday en rad m⁻².

# Retour
Trois cubes (nx, ny, nφ) : |F|, Re F, Im F, en K. Les angles de F sont
référencés à λ₀² (dérotation incomplète) — voir `lambda0_squared`.

# Voir aussi
[`getRMSF`](@ref), [`rmsynthesis_parameters`](@ref).

# Exemple
```jldoctest
julia> absF, _, _ = RMSynthesis(Q, U, range(1.0e9, 2.0e9, 64), -100:1.0:100);
```
"""
```

Points systématiques : signature, physique + référence biblio avec équations, **unités de chaque argument**, valeur de retour, invariants (croissance, dimensions), `# Voir aussi`, exemple `jldoctest` (testé par CI).

### 8.3 Manuel & tutoriel (Documenter.jl)

```
docs/
├── index.md            Qu'est-ce que MOOSE, installation, citation
├── tutorial/
│   ├── 01_quickstart.md     première observation en 10 lignes (cube jouet fourni)
│   ├── 02_config.md         le JSON champ par champ + presets
│   ├── 03_faraday.md        RM, synthèse RM, lire une FDF, RMSF
│   └── 04_cluster.md        CLI, SLURM, dry-run, mémoire
├── manual/
│   ├── conventions.md       ⭐ repères LOS, signe du RM, angle IAU, unités
│   ├── physics.md           émissivité Padovani+21, Wolfire+03, hypothèses
│   ├── outputs.md           chaque fichier FITS : contenu, unités, en-têtes
│   └── performance.md       threads, Float32, tailles de cubes
└── api/                  docstrings auto (Documenter @autodocs par sous-module)
```

La page **conventions.md est la plus importante du projet** : chaque bug critique trouvé (BUG-1, 3, 6, 13) est une convention non écrite.

---

## 9. Priorisation

### 🔴 Critique — résultats scientifiques faux (à corriger avant TOUTE diffusion)

| # | Problème | Fichier | Effet |
|---|---|---|---|
| 1 | Mapping B1/B2 anti-cyclique en LOS y (BUG-1) | `ReadSimulation.jl` | Signe de U faux, angles incohérents entre LOS |
| 2 | SNR utilisé comme σ du bruit (BUG-2) | `ProcessSynchrotron.jl` | Niveau de bruit sans rapport avec la demande |
| 3 | `fftfreq(n, Δx)` au lieu de `1/Δx` + unités px/pc (BUG-3) | `Filter.jl` | Échelles filtrées fausses |
| 4 | Reshape transposé de la table d'émissivité (BUG-4) | `QUnu.jl`, `Tnu.jl` | Émissivités permutées si table ≠ carrée/ordre attendu |
| 5 | Course de threads sur le cache d'en-têtes (BUG-5) | `WriteDataOnDisk.jl` | Crashs/headers corrompus intermittents |
| 6 | Unités d'en-têtes FITS fausses (BUG-6, MHz/"Hz", BUNITs) | `DictHeaderParameters.jl` | Fichiers publics inutilisables par des tiers |

Chacun doit être accompagné d'un **test de régression** (LOS croisées, sinus filtré, table rectangulaire, relecture des en-têtes) — sinon ils reviendront.

### 🟠 Important — robustesse et performance majeures

7. Arguments fantômes + `BoxLength_pix` sans effet, validation taille cube/config (BUG-7).
8. `Tnu3D` : interpolateur par pixel/fréquence (BUG-8) — ×50+ de gain.
9. `read_FITS_file` sans validation ni fermeture (BUG-9).
10. Gardes NaN/Inf : PolarizationFraction, Borientation, Wolfire (BUG-10).
11. `get(emiss_cfg, "path")` (BUG-11) et unification Hz internes (BUG-12).
12. Brancher `rmsynthesis_parameters`/`getRMSF` : écrire RMSF + FWHM + φ_max (§3.1).
13. Pipeline Float32 + suppression des temporaires (§4.1) — condition d'accès au 512³ confortable.
14. Validation précoce + estimation mémoire + `--dry-run` (§6.2, §7.6).

### 🟡 Moyen — qualité, maintenabilité, science fine

15. Convention cumsum RM documentée + point milieu (BUG-13).
16. Interpolation d'émissivité en log-log + politique d'extrapolation journalisée (§3.3).
17. Refactor `NeModel`/RunConfig typée/fusion des 3 ProcessSynchrotron (§5.3).
18. Poids dans RMSynthesis ; `LAMBDA0SQ` dans l'en-tête (§3.1).
19. Apodisation + normalisation documentée des spectres de puissance (§3.6).
20. ProgressMeter avec ETA, Logging structuré, résumé enrichi (§7).
21. Docstrings du cœur scientifique + Documenter.jl + conventions.md (§8).
22. Suppression code mort & dépendances fantômes ; `Comonicon`/`ArgParse` pour le CLI.

### ⚪ Faible

23. Nommage snake_case + `@deprecate` ; logo conditionné à `isatty` ; étiquettes de tracés en anglais ; `mv` T_nu ; commentaire de `C_m` ; double séparateur CreateFreqFile ; mots-clés FITS DATE/ORIGIN/VERSION.

**Séquencement recommandé** : (i) écrire d'abord les tests de régression qui *capturent le comportement actuel correct* (LOS z, RM, Tb, synthèse RM sur cas analytiques) ; (ii) corriger les 6 critiques avec leurs tests ; (iii) publier une version `v0.x` corrigée avec CHANGELOG explicitant les changements de résultats (LOS y, bruit, filtre) ; (iv) dérouler Important puis Moyen.

---

## 10. Vision long terme

### 10.1 Fonctionnalités manquantes utiles à la communauté

- **RMSF systématique + RM-CLEAN** (déconvolution de la FDF) — sans cela les FDF restent qualitatives ; interop avec RM-Tools (sorties compatibles).
- **Dépolarisation intra-cellule** (facteur de Burn) en option (§3.2).
- **Poids/canaux drapeautés** dans la synthèse RM (RFI réalistes).
- **Convolution par un lobe** (beam gaussien) et grilles de visibilités simplifiées — l'actuel masque 0/1 est une approximation grossière d'un interféromètre ; à terme, export vers des simulateurs (OSKAR, pyuvsim).
- **Mode tranches / out-of-core** pour 2048³ (§4.1.5) et, au-delà, distribution multi-nœuds (Distributed.jl ou MPI.jl) — la séparabilité en colonnes rend cela naturel.
- **Unités explicites** : Unitful.jl aux frontières de l'API publique (entrées converties puis strippées en interne pour la performance) — élimine structurellement la classe de bugs MHz/Hz/pc/cm.
- Sorties **HDF5/Zarr** optionnelles pour les très gros cubes (FITS reste le format d'échange).

### 10.2 Standards & paquets de l'écosystème à adopter

| Besoin | Paquet |
|---|---|
| CLI | Comonicon.jl ou ArgParse.jl |
| Progression/ETA | ProgressMeter.jl |
| Config validée | Configurations.jl ou JSON3 + StructTypes |
| Analyse statique | JET.jl, Aqua.jl (déps fantômes, ambiguïtés, piracy) |
| Bench | BenchmarkTools.jl + PkgBenchmark.jl |
| Doc | Documenter.jl + doctests en CI |
| Unités | Unitful.jl, UnitfulAstro.jl |
| WCS/headers | WCS.jl (en-têtes FITS standards vérifiés) |

### 10.3 Pratiques de développement

- **CI GitHub Actions** : tests sur Julia LTS + stable, Linux/macOS, 1 et 4 threads (BUG-5 n'aurait pas survécu à une CI multi-thread) ; Aqua + JET + doctests ; couverture Codecov.
- **SemVer + CHANGELOG** ; enregistrement au **General registry** Julia ; `CITATION.cff` (la citation 2026A&A...708A.245B existe déjà dans le README — la rendre machine-lisible) ; licence claire ; CONTRIBUTING.md.
- **Revue par PR systématique**, même seul (auto-revue à froid) ; branches de fonctionnalités ; jamais de push direct sur main.

### 10.4 Tests unitaires et de validation scientifique (le manque le plus grave après les bugs)

Le `runtests.jl` actuel ne teste **aucune physique**. Suite cible :

1. **Analytiques exacts** : champ B uniforme + ne uniforme → RM linéaire en z connu ; ψ_src constant → Q/U = ε_P·n_z·(cos, sin)(2χ) vérifiable à la machine près ; Tb d'une émissivité constante.
2. **Synthèse RM** : écran de Faraday mince à φ₀ injecté analytiquement → pic de |F| à φ₀, FWHM ≈ 3.8/Δλ², et `RMSynthesis(RMSF du signal unité) == getRMSF`.
3. **Invariance LOS** : cube isotrope statistiquement → moments de I identiques sur x/y/z ; dipôle ponctuel → position et signes attendus (verrouille BUG-1 à jamais).
4. **Filtre** : sinus pur conservé/annihilé selon la bande (verrouille BUG-3).
5. **Émissivité** : table rectangulaire synthétique f(B,ν) connue → interpolateur exact aux nœuds, erreur bornée entre nœuds (verrouille BUG-4) ; comparaison ponctuelle interpolation vs intégration QuadGK directe (< 1 % en log-log).
6. **Bruit** : σ mesuré a posteriori ≈ signal/SNR demandé (verrouille BUG-2) ; reproductibilité par graine.
7. **Round-trip FITS** : chaque produit relu → mots-clés d'unités conformes (verrouille BUG-6) ; lecture par astropy dans un test d'intégration optionnel.
8. **Propriétés** (à terme, Supposition.jl) : RM(cumsum) croissant si ne·B∥ > 0 ; |p| ≤ 1 ; F(φ) hermitienne pour U = 0.

### 10.5 Cap proposé

À 12 mois : MOOSE enregistré, CI verte multi-thread, doc Documenter publiée avec tutoriel exécutable, suite de validation physique, pipeline Float32 streaming gérant 1024³ sur une station 64 Gio — et un CHANGELOG honnête documentant la correction des résultats LOS y/bruit/filtre. C'est la différence entre « un code de groupe » et « un outil communautaire citable ».

---

## Synthèse exécutive

MOOSE a un **cœur scientifique majoritairement sain** (synthèse RM conforme à B&dB 2005, prefacteur RM correct, décomposition Q/U standard, Tb correcte) et une vraie richesse fonctionnelle (4 modèles de ne, HEALPix, filtrage, statistiques). Mais la diffusion publique est aujourd'hui **bloquée par six défauts critiques** qui produisent des résultats faux : le repère LOS y inversé, le bruit dont l'amplitude croît avec le SNR demandé, le filtre aux échelles fausses, la table d'émissivité transposée, une course de threads à l'écriture, et des en-têtes FITS aux unités erronées. Tous sont corrigeables en quelques jours ; aucun correctif ne supprime de fonctionnalité ; trois d'entre eux (LOS y, bruit, filtre) **changent les résultats — parce que les résultats actuels sont erronés** — et doivent être documentés comme tels. La priorité absolue n'est pas d'optimiser ni de refactorer, mais d'écrire les tests de validation physique qui n'existent pas, puis de corriger sous leur protection.

