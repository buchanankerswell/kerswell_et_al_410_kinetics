#######################################################
## .0. Load Libraries                            !!! ##
#######################################################
from pathlib import Path

from burnman import minerals
from profiles import AdiabaticProfile, DrivingForceProfile, parse_arguments


#######################################################
## .1. Main                                      !!! ##
#######################################################
def main() -> None:
    """Main function to generate ASPECT profiles."""
    args = parse_arguments()

    out_dir = Path(args.out_dir)
    out_fig_dir = Path(args.out_fig_dir)
    out_table = out_dir / f"{args.model_id}-material-table.txt"
    out_adiabat = out_dir / f"{args.model_id}-adiabatic-profile.txt"
    out_driving_force = out_dir / "thermodynamic-driving-force-profile-olivine-wadsleyite.txt"

    adiabat = AdiabaticProfile(
        model_id=args.model_id,
        out_table=out_table,
        out_profile=out_adiabat,
        potential_temperature=args.potential_temperature,
        out_table_resolution=args.out_table_resolution,
        out_profile_resolution=args.out_profile_resolution,
        planet_radius=args.planet_radius,
        surface_gravity=args.surface_gravity,
    )

    ol = minerals.SLB_2024.olivine(molar_fractions=[0.8, 0.2])
    wad = minerals.SLB_2024.wadsleyite(molar_fractions=[0.8, 0.2])

    driving_force = DrivingForceProfile(
        material_a=ol,
        material_b=wad,
        in_profile=adiabat,
        out_profile=out_driving_force,
        out_fig_dir=out_fig_dir,
    )
    driving_force.visualize()
    driving_force._visualize_delta_PT_at_target_depth(target_depth=420e3)
    driving_force._visualize_delta_PT_at_target_depth(target_depth=430e3)
    driving_force._visualize_delta_PT_at_target_depth(target_depth=440e3)
    driving_force._visualize_delta_PT_at_target_depth(target_depth=450e3)


if __name__ == "__main__":
    main()
