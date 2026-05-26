!-------------------------------------------------------------------------------
! EXTERNAL NUTRIENT SINK FOR KELP FARM SIMULATION (DKS)
!-------------------------------------------------------------------------------
!
! This implementation provides a one-way coupled kelp farm simulation where
! nutrient and carbon uptake rates are prescribed from a standalone kelp growth
! model and applied as sinks to the COBALT biogeochemical model.
!
! APPROACH:
! Uptake rates are computed offline by a standalone kelp growth model and
! prescribed to COBALT via the FMS data_override mechanism. The kelp biomass
! is assumed to be harvested, so removed nutrients and carbon do not re-enter
! the ocean through detrital or dissolved organic pools.
!
! LIMITATIONS:
! This is a one-way coupling - the kelp growth model does not see the nutrient
! depletion caused by the farm in COBALT. As the simulation progresses, local
! nutrient concentrations in COBALT may diverge from those assumed by the kelp
! model. A hard limiter prevents removal of more nutrients than are available
! in a given grid cell, but uptake rates are not gracefully reduced as nutrients
! decline. This approach is most appropriate for:
!   - Short simulations
!   - Low-density farms
!   - Nutrient-rich environments where depletion is minimal
!
! REQUIRED INPUTS (via data_override, field name in ocean model: 'OCN'):
!
!   mask_e_juptake (dimensionless, 3D: isc:iec, jsc:jec, 1:nk)
!     - Spatial and vertical mask defining the extent of the kelp farm
!     - Values > 0 activate the sink at that location
!     - Should be restricted to the photic zone and farm footprint
!     - Units: dimensionless
!
!   e_juptake_no3 (mol NO3 kg-1 s-1, 3D: isc:iec, jsc:jec, 1:nk)
!     - Nitrate uptake rate prescribed from the kelp growth model
!     - Primary nitrogen source for kelp in the model
!     - NH4 uptake is neglected (assumed small relative to NO3)
!     - Applied with limiter: cannot exceed available NO3 concentration
!
!   e_juptake_po4 (mol PO4 kg-1 s-1, 3D: isc:iec, jsc:jec, 1:nk)
!     - Phosphate uptake rate prescribed from the kelp growth model
!     - Should be stoichiometrically consistent with e_juptake_no3
!     - Applied with limiter: cannot exceed available PO4 concentration
!
!   e_juptake_fed (mol Fe kg-1 s-1, 3D: isc:iec, jsc:jec, 1:nk)
!     - Dissolved iron uptake rate prescribed from the kelp growth model
!     - Applied with limiter: cannot exceed available Fed concentration
!
! STOICHIOMETRIC ADJUSTMENTS APPLIED INTERNALLY:
!
!   NO3:  decremented by e_juptake_no3
!   PO4:  decremented by e_juptake_po4
!   Fed:  decremented by e_juptake_fed
!   DIC:  decremented by c_2_n * e_juptake_no3  (Redfield C:N = 106/16)
!   O2:   incremented by o2_2_no3 * e_juptake_no3 (= 150/16 mol O2 mol N-1)
!   ALK:  incremented by e_juptake_no3  (NO3 uptake increases alkalinity
!                                        by 1 equivalent per mole NO3)
!
! NOTE ON ALK: If NH4 uptake were included it would decrease alkalinity,
! partially offsetting the NO3-driven increase. Neglecting NH4 uptake
! therefore leads to a small overestimate of the alkalinity increase,
! which is considered acceptable given the assumed small NH4 contribution.
!
! MASS BALANCE:
! The external uptake is accounted for in the pre/post source-sink mass
! balance checker by subtracting e_juptake_no3, e_juptake_po4 and
! e_juptake_fed from pre_totn, pre_totp and pre_totfe respectively,
! and subtracting c_2_n*e_juptake_no3 from pre_totc before the
! imbalance check is performed.
!
!--------------------------------------------------------------------------------```