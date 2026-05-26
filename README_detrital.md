!-------------------------------------------------------------------------
! EXTERNAL ORGANIC MATTER ADDITION TO SEAFLOOR (DKS)
!-------------------------------------------------------------------------
!
! This implementation prescribes an addition of particulate organic matter
! directly to the bottom layer of the ocean at specified locations. It is
! designed to be used independently of the external nutrient sink block
! (Section 4.5). If both blocks are run simultaneously the user is responsible
! for ensuring that the prescribed fluxes are stoichiometrically consistent
! to avoid spurious creation or destruction of mass at the farm scale.
!
! APPROACH:
! Organic matter addition rates are prescribed via the FMS data_override
! mechanism and injected into the detrital pools (ndet, pdet, fedet) at
! the bottom layer of each active grid cell. The surface level of the 3D
! mask (mask_addition_t(i,j,1)) is used as a 2D horizontal on/off switch
! to identify the farm footprint, and material is deposited at
! k = grid_kmt(i,j) regardless of depth. This simulates instantaneous
! sinking of organic material to the seafloor (e.g., harvested kelp
! biomass or natural senescent kelp detritus).
!
! DIC and ALK are not adjusted at the moment of injection. The carbon
! cycle impact occurs downstream when the injected detritus remineralizes
! via COBALT's existing remineralization machinery (jremin_ndet), which
! handles DIC and ALK adjustments automatically.
!
! LIMITATIONS:
! - Sinking is instantaneous: water column remineralization during descent
!   is bypassed entirely. This is most appropriate for shallow farm
!   locations or rapidly sinking material.
! - There is no limiter on the addition rate. Large prescribed fluxes
!   relative to the bottom layer thickness can produce unrealistically
!   large detrital concentrations in a single timestep. The user should
!   verify that prescribed fluxes are physically reasonable relative to
!   local grid cell volumes.
! - The three detrital pools (ndet, pdet, fedet) are prescribed
!   independently. Stoichiometric consistency between them is the
!   responsibility of the user.
!
! REQUIRED INPUTS (via data_override, field name in ocean model: 'OCN'):
!
!   mask_addition_t (dimensionless, 3D: isc:iec, jsc:jec, 1:nk)
!     - Only the surface level mask_addition_t(i,j,1) is used
!     - Values > 0 activate the addition at that horizontal location
!     - Material is always deposited at k = grid_kmt(i,j)
!     - Units: dimensionless
!
!   ndet_addition (mol N kg-1 s-1, 2D: isc:iec, jsc:jec)
!     - Detrital nitrogen addition rate
!     - Applied to p_ndet at the bottom layer
!
!   pdet_addition (mol P kg-1 s-1, 2D: isc:iec, jsc:jec)
!     - Detrital phosphorus addition rate
!     - Applied to p_pdet at the bottom layer
!     - Should be stoichiometrically consistent with ndet_addition
!
!   fedet_addition (mol Fe kg-1 s-1, 2D: isc:iec, jsc:jec)
!     - Detrital iron addition rate
!     - Applied to p_fedet at the bottom layer
!     - Should be stoichiometrically consistent with ndet_addition
!
! TRACER ADJUSTMENTS APPLIED INTERNALLY:
!
!   ndet:  incremented by ndet_addition
!   pdet:  incremented by pdet_addition
!   fedet: incremented by fedet_addition
!   DIC:   not adjusted at injection (handled by remineralization)
!   ALK:   not adjusted at injection (handled by remineralization)
!
! MASS BALANCE:
! The additions are accounted for in the pre/post source-sink mass balance
! checker by incrementing pre_totn, pre_totp, pre_totfe and pre_totc
! (the latter via c_2_n * ndet_addition) before the imbalance check
! is performed.
!
!
!-------------------------------------------------------------------------