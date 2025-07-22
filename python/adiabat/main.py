#######################################################
## .0. Load Libraries                            !!! ##
#######################################################
from pathlib import Path

from burnman import minerals, Composite
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
    out_forward_reaction_410 = out_dir / "forward-reaction-410-profile.txt"
    out_reverse_reaction_410 = out_dir / "reverse-reaction-410-profile.txt"
    out_forward_reaction_550 = out_dir / "forward-reaction-550-profile.txt"
    out_reverse_reaction_550 = out_dir / "reverse-reaction-550-profile.txt"
    out_forward_reaction_660 = out_dir / "forward-reaction-660-profile.txt"
    out_reverse_reaction_660 = out_dir / "reverse-reaction-660-profile.txt"

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

    # Define minerals
    ol = minerals.SLB_2024.olivine([1.0, 0.0])
    wa = minerals.SLB_2024.wadsleyite([1.0, 0.0])
    ri = minerals.SLB_2024.ringwoodite([1.0, 0.0])
    pv = minerals.SLB_2024.bridgmanite([1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    fper = minerals.SLB_2024.ferropericlase([1.0, 0.0, 0.0, 0.0, 0.0])

    # Define mantle materials
    upper_mantle = Composite([ol], [1.0], name="upper-mantle")
    upper_mtz = Composite([wa], [1.0], name="upper-mtz")
    lower_mtz = Composite([ri], [1.0], name="lower-mtz")
    lower_mantle = Composite([pv, fper], [0.8, 0.2], name="lower-mantle")

    forward_reaction_410 = DrivingForceProfile(
        material_a=upper_mantle,
        material_b=upper_mtz,
        in_profile=adiabat,
        out_profile=out_forward_reaction_410,
        out_fig_dir=out_fig_dir,
    )
    forward_reaction_410.visualize()

    reverse_reaction_410 = DrivingForceProfile(
        material_a=upper_mtz,
        material_b=upper_mantle,
        in_profile=adiabat,
        out_profile=out_reverse_reaction_410,
        out_fig_dir=out_fig_dir,
    )
    reverse_reaction_410.visualize()

    forward_reaction_550 = DrivingForceProfile(
        material_a=upper_mtz,
        material_b=lower_mtz,
        in_profile=adiabat,
        out_profile=out_forward_reaction_550,
        out_fig_dir=out_fig_dir,
    )
    forward_reaction_550.visualize()

    reverse_reaction_550 = DrivingForceProfile(
        material_a=lower_mtz,
        material_b=upper_mtz,
        in_profile=adiabat,
        out_profile=out_reverse_reaction_550,
        out_fig_dir=out_fig_dir,
    )
    reverse_reaction_550.visualize()

    forward_reaction_660 = DrivingForceProfile(
        material_a=lower_mtz,
        material_b=lower_mantle,
        in_profile=adiabat,
        out_profile=out_forward_reaction_660,
        out_fig_dir=out_fig_dir,
    )
    forward_reaction_660.visualize()

    reverse_reaction_660 = DrivingForceProfile(
        material_a=lower_mantle,
        material_b=lower_mtz,
        in_profile=adiabat,
        out_profile=out_reverse_reaction_660,
        out_fig_dir=out_fig_dir,
    )
    reverse_reaction_660.visualize()


if __name__ == "__main__":
    main()
