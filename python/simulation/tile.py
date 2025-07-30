#######################################################
## .0. Load Libraries                            !!! ##
#######################################################
import platform
import re
from dataclasses import dataclass, field
from pathlib import Path

from moviepy.video.io.ImageSequenceClip import ImageSequenceClip
from PIL import Image, ImageDraw, ImageFont
from visualize import PyVistaModelConfig


#######################################################
## .1. ImageTiler                                !!! ##
#######################################################
@dataclass
class ImageTiler:
    plot_config: PyVistaModelConfig
    out_fig_dir: Path
    field1: str
    field2: str
    field3: str | None = None
    movie: bool = False
    fps: int = 10
    tags: str | None = None
    tag_size: int = 188

    tiles_dir: Path = field(init=False)
    movies_dir: Path = field(init=False)
    prefix: str = field(init=False)
    images_field1: list[Path] = field(init=False)
    images_field2: list[Path] = field(init=False)
    images_field3: list[Path] = field(init=False)
    tiled_images: list[str] = field(default_factory=list, init=False)

    verbosity: int = 0

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __post_init__(self):
        """"""
        self.tiles_dir = self.out_fig_dir / "tiles"
        self.movies_dir = self.out_fig_dir / "movies"

        self.prefix = str(self.out_fig_dir.name).replace("_", "-")
        self.images_field1 = self._get_images(self.field1)
        self.images_field2 = self._get_images(self.field2)
        self.images_field3 = self._get_images(self.field3) if self.field3 else []

        fields_list = [self.field1, self.field2]
        if self.field3:
            fields_list.append(self.field3)

        self.field_file_mappings = [
            self.plot_config.file_mapping.get(field, None) for field in fields_list
        ]

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def _detect_prefix(self) -> str:
        """"""
        sample = next(self.out_fig_dir.glob(f"*{self.field1}*.png"), None)
        if sample:
            return re.sub(rf"{self.field1}-.*\.png", "", sample.name).rstrip("-")
        raise ValueError(
            f"No matching images found for variable '{self.field1}' to detect prefix."
        )

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def _get_images(self, field: str | None) -> list[Path]:
        """"""
        if not field:
            return []

        field_file_mapping = self.plot_config.file_mapping.get(field, None)

        if field_file_mapping:
            return sorted(
                self.out_fig_dir.glob(f"{self.prefix}-{field_file_mapping}-*.png")
            )
        else:
            return []

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def _get_binned_depth_images(self) -> list[Path]:
        """"""
        if not self.out_fig_dir.exists():
            return []

        all_file_mappings_none = all(
            mapping is None for mapping in self.field_file_mappings
        )

        if not all_file_mappings_none:
            pattern = (
                f"*binned*-{self.field_file_mappings[0]}-{self.field_file_mappings[1]}"
            )
            if self.field3:
                pattern += f"-{self.field_file_mappings[2]}"

            pattern += "-*.png"

            return sorted(self.out_fig_dir.glob(pattern))
        else:
            return []

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def _extract_timestep(self, filename: Path) -> str | None:
        """"""
        match = re.search(r"-(\d{4})\.png$", filename.name)
        return match.group(1) if match else None

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def _get_default_font(self) -> ImageFont.ImageFont | ImageFont.FreeTypeFont:
        """"""

        system = platform.system()

        font_paths = {
            "Darwin": "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
            "Windows": "C:/Windows/Fonts/arialbd.ttf",
            "Linux": "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
        }

        font_path = font_paths.get(system)

        if font_path and Path(font_path).exists():
            try:
                return ImageFont.truetype(font_path, self.tag_size)
            except IOError:
                if self.verbosity >= 1:
                    print(f" !! Warning: failed to load truetype font at {font_path}!")

        if self.verbosity >= 1:
            print(" !! Warning: falling back to default font (tiny and non-scalable)!")

        return ImageFont.load_default()

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def _annotate_image(self, img_path: Path, tag: str | None) -> Image.Image:
        """"""
        image = Image.open(img_path).convert("RGBA")

        if tag is None:
            return image

        font = self._get_default_font()
        draw = ImageDraw.Draw(image)
        draw.text((25, 10), tag, font=font, fill="black")

        return image

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def _concatenate_images_horizontally(
        self, images: list[Image.Image]
    ) -> Image.Image:
        """"""
        widths, heights = zip(*(img.size for img in images))
        total_width = sum(widths)
        max_height = max(heights)

        combined = Image.new("RGB", (total_width, max_height), color=(255, 255, 255))

        x_offset = 0
        for img in images:
            combined.paste(img, (x_offset, 0))
            x_offset += img.size[0]

        return combined

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def tile_images(self) -> None:
        """"""
        fields_list = [self.field1, self.field2]
        if self.field3:
            fields_list.append(self.field3)

        for i, (img1, img2) in enumerate(zip(self.images_field1, self.images_field2)):
            ts1 = self._extract_timestep(img1)
            ts2 = self._extract_timestep(img2)

            img3 = (
                self.images_field3[i]
                if self.field3 and i < len(self.images_field3)
                else None
            )

            ts3 = self._extract_timestep(img3) if img3 else ts1

            if ts1 != ts2 or ts1 != ts3:
                if self.verbosity >= 1:
                    print(" !! Warning: timestep mismatch!\n" " -- Skipping tile")
                continue

            tags = (
                [c.lower() for c in self.tags.ljust(3)[:3]]
                if self.tags is not None
                else [None, None, None]
            )

            annotated_imgs = [
                self._annotate_image(img1, tags[0]),
                self._annotate_image(img2, tags[1]),
            ]

            if img3:
                annotated_imgs.append(self._annotate_image(img3, tags[2]))

            all_tags_none = all(tag is None for tag in tags)

            out_tile = (
                self.tiles_dir
                / f"{self.prefix}-{self.field_file_mappings[0]}-{self.field_file_mappings[1]}"
                f"{'-' + str(self.field_file_mappings[2]) if self.field_file_mappings[2] is not None else ''}"
                f"-{ts1}"
                f"{'-tagged' if not all_tags_none else ''}.png"
            )

            if out_tile.exists():
                print(f" -- Found tile: {out_tile.name}!")
                continue

            self.tiles_dir.mkdir(parents=True, exist_ok=True)
            combined = self._concatenate_images_horizontally(annotated_imgs)
            combined.save(out_tile)
            combined.close()

            self.tiled_images.append(str(out_tile))

            print(f"--> Tile: {out_tile.name}")

        if self.movie and self.tiled_images:
            self._create_movie()

        self._compose_with_binned_depth()

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def _compose_with_binned_depth(self) -> None:
        """"""
        binned_images = self._get_binned_depth_images()
        if not binned_images:
            if self.verbosity >= 1:
                print(" !! Warning: no binned depth images found!")
            return

        composition_dir = self.out_fig_dir / "composition"
        composition_dir.mkdir(parents=True, exist_ok=True)

        composed_tiles = []

        for tile_path, binned_path in zip(self.tiled_images, binned_images):
            ts_tile = self._extract_timestep(Path(tile_path))
            ts_binned = self._extract_timestep(Path(binned_path))

            if ts_tile != ts_binned:
                if self.verbosity >= 1:
                    print(
                        f" !! Warning: timestep {ts_tile} mismatch!\n -- Skipping tile"
                    )
                continue

            img_tile = Image.open(tile_path).convert("RGB")
            img_binned = Image.open(binned_path).convert("RGB")

            # Pad binned image to match width
            width_tile, _ = img_tile.size
            width_binned, height_binned = img_binned.size
            padded_binned = Image.new("RGB", (width_tile, height_binned), "white")
            offset = ((width_tile - width_binned) // 2, 0)
            padded_binned.paste(img_binned, offset)

            # Concatenate vertically
            composed = Image.new(
                "RGB", (width_tile, img_tile.height + padded_binned.height), "white"
            )
            composed.paste(img_tile, (0, 0))
            composed.paste(padded_binned, (0, img_tile.height))

            out_path = composition_dir / Path(tile_path).name
            composed.save(out_path)
            composed.close()

            composed_tiles.append(str(out_path))
            print(f"--> Composition: {out_path.name}")

        if self.movie and composed_tiles:
            movie_path = (
                self.out_fig_dir
                / "movies_comps"
                / f"{self.prefix}-{self.field_file_mappings[0]}-{self.field_file_mappings[1]}{'-' + self.field_file_mappings[2] if self.field_file_mappings[2] else ''}.mp4"
            )
            movie_path.parent.mkdir(exist_ok=True, parents=True)

            clip = ImageSequenceClip(composed_tiles, fps=self.fps)
            clip.write_videofile(
                str(movie_path),
                codec="libx264",
                audio=False,
                verbose=False,
                logger=None,
            )
            clip.close()

            print(f"--> Composition movie: {movie_path.name}")

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def _create_movie(self) -> None:
        """"""
        out_path = (
            self.movies_dir
            / f"{self.prefix}-{self.field_file_mappings[0]}-{self.field_file_mappings[1]}{'-' + self.field_file_mappings[2] if self.field_file_mappings[2] else ''}.mp4"
        )

        if out_path.exists():
            print(f" -- Found movie: {out_path.name}!")
            return

        clip = ImageSequenceClip(self.tiled_images, fps=self.fps)
        clip.write_videofile(
            str(out_path), codec="libx264", audio=False, verbose=False, logger=None
        )
        clip.close()

        print(f"--> Movie: {out_path.name}")
