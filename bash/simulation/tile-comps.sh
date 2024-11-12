#!/usr/bin/env bash

if [ $# -lt 3 ]; then
    echo "Usage: $0 <image_dir> <variable_name_1> <variable_name_2> [variable_name_3] [--movie=<true|false>] [--fps=<fps>] [--tag-char=<abc>]"
    exit 1
fi

IMG_DIR=$1
VAR1=$2
VAR2=$3
VAR3=""

MOVIE="false"
FPS=10
TAG_CHAR="abc"
TAG_SIZE=164

shift 3
while (( "$#" )); do
    case "$1" in
        --movie=*)
            MOVIE="${1#*=}"
            ;;
        -movie)
            shift
            MOVIE=$1
            ;;
        -m)
            shift
            MOVIE=$1
            ;;
        --fps=*)
            FPS="${1#*=}"
            ;;
        -fps)
            shift
            FPS=$1
            ;;
        -f)
            shift
            FPS=$1
            ;;
        --tag-char=*)
            TAG_CHAR="${1#*=}"
            ;;
        -tag-char)
            shift
            TAG_CHAR=$1
            ;;
        -t)
            shift
            TAG_CHAR=$1
            ;;
        *)
            VAR3=$1
            ;;
    esac
    shift
done

if [[ "$MOVIE" != "true" && "$MOVIE" != "false" ]]; then
    echo "Error: Movie must be either 'true' or 'false' ..."
    exit 1
fi

if ! [[ "$FPS" =~ ^[0-9]+$ ]]; then
    echo "Error: FPS must be a positive integer ..."
    exit 1
fi

SAMPLE=$(find "$IMG_DIR" -name "smp-pda-${VAR1}*.png" | head -n 1)
PREFIX=$(echo "$SAMPLE" | sed -r "s|^${IMG_DIR}/(.*)-${VAR1}-.*\.png|\1|")

mapfile -t IMGS_VAR1 < <(find "$IMG_DIR" -name "${PREFIX}-${VAR1}-*.png" | sort -u)
mapfile -t IMGS_VAR2 < <(find "$IMG_DIR" -name "${PREFIX}-${VAR2}-*.png" | sort -u)

IMGS_VAR3=()
if [ -n "$VAR3" ]; then
    mapfile -t IMGS_VAR3 < <(find "$IMG_DIR" -name "${PREFIX}-${VAR3}-*.png" | sort -u)
fi

TILED_IMAGES=()

echo "========================================================="
echo "Tiling compositions with imagemagick ..."
echo "========================================================="
echo "Image directory: $IMG_DIR"
echo "Tiles directory: $IMG_DIR/tiles_comps"
echo "Movie directory: $IMG_DIR/movies_comps"
echo "Variable 1:      $VAR1"
echo "Variable 2:      $VAR2"
echo "Variable 3:      $VAR3"
echo "Movie:           $MOVIE"
echo "N images:        ${#IMGS_VAR1[@]}"
echo "Framerate:       $FPS"
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "Tiling begins ..."
echo "---------------------------------------------------------"

for ((i = 1; i <= ${#IMGS_VAR1[@]}; i++)); do
    TILE_COUNTER=1

    IMG1="${IMGS_VAR1[$i]}"
    IMG2="${IMGS_VAR2[$i]}"
    IMG3=""
    if [[ -n "$VAR3" ]]; then
        IMG3="${IMGS_VAR3[$i]}"
    fi
    T1=$(basename "$IMG1" ".png" | sed -r "s/^.*-([0-9]{4})$/\1/")
    T2=$(basename "$IMG2" ".png" | sed -r "s/^.*-([0-9]{4})$/\1/")
    T3=$T1
    if [[ -n "$IMG3" ]]; then
        T3=$(basename "$IMG3" ".png" | sed -r "s/^.*-([0-9]{4})$/\1/")
    fi

    mkdir -p "${IMG_DIR}/tiles"
    OUTPUT_TILE="${IMG_DIR}/tiles/${VAR1}-${VAR2}${VAR3:+-$VAR3}-${T1}.png"

    if [[ -f "$OUTPUT_TILE" ]]; then
        TILED_IMAGES+=("$OUTPUT_TILE")
    elif [[ -f "$IMG1" && -f "$IMG2" && "$T1" == "$T2" && "$T1" == "$T3" ]]; then
        echo -ne "Processing images for timestep $T1 ... \r"

        IMG1_TAG="${IMG_DIR}/${VAR1}-${VAR2}${VAR3:+-$VAR3}-${T1}-tag1.png"
        IMG2_TAG="${IMG_DIR}/${VAR1}-${VAR2}${VAR3:+-$VAR3}-${T1}-tag2.png"
        IMG3_TAG="${IMG_DIR}/${VAR1}-${VAR2}${VAR3:+-$VAR3}-${T1}-tag3.png"

        TAG="$(echo "$TILE_COUNTER" | awk -v tag_char="$TAG_CHAR" \
            '{print tolower(substr(tag_char, $1, 1))}')"
        magick "$IMG1" -gravity NorthWest -pointsize "$TAG_SIZE" \
            -annotate +15+10 "$TAG" "$IMG1_TAG"
        ((TILE_COUNTER++))

        TAG="$(echo "$TILE_COUNTER" | awk -v tag_char="$TAG_CHAR" \
            '{print tolower(substr(tag_char, $1, 1))}')"
        magick "$IMG2" -gravity NorthWest -pointsize "$TAG_SIZE" \
            -annotate +15+10 "$TAG" "$IMG2_TAG"
        ((TILE_COUNTER++))

        if [[ -n "$IMG3" ]]; then
            TAG="$(echo "$TILE_COUNTER" | awk -v tag_char="$TAG_CHAR" \
                '{print tolower(substr(tag_char, $1, 1))}')"
            magick "$IMG3" -gravity NorthWest -pointsize "$TAG_SIZE" \
                -annotate +15+10 "$TAG" "$IMG3_TAG"
            ((TILE_COUNTER++))
            magick "$IMG1_TAG" "$IMG2_TAG" "$IMG3_TAG" \
                +append "$OUTPUT_TILE"
        else
            magick "$IMG1_TAG" "$IMG2_TAG" +append "$OUTPUT_TILE"
        fi

        rm "$IMG1_TAG" "$IMG2_TAG" "$IMG3_TAG"
#        pngquant --quality=80-95 --speed=1 --floyd=0 --ext .png --force "$OUTPUT_TILE"

        TILED_IMAGES+=("$OUTPUT_TILE")
    else
        echo "Skipping timestep $T1: Images don't exist or tsteps don't match."
    fi
done

mapfile -t BINNED_IMAGES < <(find "$IMG_DIR/binned_depth" -name "*${VAR1}-${VAR2}${VAR3:+-$VAR3}*.png" | sort -u)

if [ ${#BINNED_IMAGES[@]} -eq 0 ]; then
    echo "No binned images found!"
    exit 1
fi

TILED_COMPS=()

for ((i = 1; i <= ${#TILED_IMAGES[@]}; i++)); do
    TILE_COUNTER=1

    TILE1="${TILED_IMAGES[$i]}"
    TILE2="${BINNED_IMAGES[$i]}"
    T1=$(basename "$TILE1" ".png" | sed -r "s/^.*-([0-9]{4})$/\1/")
    T2=$(basename "$TILE2" ".png" | sed -r "s/^.*-([0-9]{4})$/\1/")

    mkdir -p "${IMG_DIR}/tiles_comps"
    OUTPUT_TILE="${IMG_DIR}/tiles_comps/${VAR1}-${VAR2}${VAR3:+-$VAR3}-${T1}.png"

    if [[ -f "$OUTPUT_TILE" ]]; then
        TILED_COMPS+=("$OUTPUT_TILE")
    elif [[ -f "$TILE1" && -f "$TILE2" && "$T1" == "$T2" ]]; then
        echo -ne "Processing composition for timestep $T1 ... \r"

        WIDTH_TOP=$(magick identify -format "%w" "$TILE1")
        HEIGHT_BOTTOM=$(magick identify -format "%h" "$TILE2")

        TILE2_PADDED="${IMG_DIR}/${VAR1}-${VAR2}${VAR3:+-$VAR3}-${T1}-tmp.png"

        magick "$TILE2" -background white -gravity center \
					-extent "${WIDTH_TOP}x${HEIGHT_BOTTOM}+50+0" "$TILE2_PADDED"

        magick "$TILE1" "$TILE2_PADDED" -append "$OUTPUT_TILE"

        WIDTH=$(magick identify -format "%w" "$OUTPUT_TILE")
        HEIGHT=$(magick identify -format "%h" "$OUTPUT_TILE")
        NEW_WIDTH=$(( (WIDTH + 1) & ~1 ))
        NEW_HEIGHT=$(( (HEIGHT + 1) & ~1 ))

        if [[ $WIDTH -ne $NEW_WIDTH || $HEIGHT -ne $NEW_HEIGHT ]]; then
            magick "$OUTPUT_TILE" -background white -gravity center \
                   -extent "${NEW_WIDTH}x${NEW_HEIGHT}" "$OUTPUT_TILE"
        fi

#        pngquant --quality=80-95 --speed=1 --floyd=0 --ext .png --force "$OUTPUT_TILE"

        rm "$TILE2_PADDED"

        ((TILE_COUNTER++))
        TILED_COMPS+=("$OUTPUT_TILE")
    else
        echo "Skipping timestep $T1: Images don't exist or tsteps don't match."
    fi

done

if [[ "$MOVIE" == "true" && ${#TILED_COMPS[@]} -gt 0 ]]; then
    mkdir -p "${IMG_DIR}/movies_comps"
    MOVIE_OUTPUT="${IMG_DIR}/movies_comps/${VAR1}-${VAR2}${VAR3:+"-$VAR3"}"

    SUFFIX_COUNTER=1
    ORIGINAL_MOVIE_OUTPUT="$MOVIE_OUTPUT"
    while [[ -f "$MOVIE_OUTPUT" ]]; do
        MOVIE_OUTPUT="${ORIGINAL_MOVIE_OUTPUT%.png}-$SUFFIX_COUNTER.png"
        SUFFIX_COUNTER=$((SUFFIX_COUNTER + 1))
    done

    if [ -f "$MOVIE_OUTPUT.mp4" ]; then
        echo "========================================================="
        echo "Movie found!"
    else
        echo "---------------------------------------------------------"
        echo "Creating movie with ffmpeg ..."
        ffmpeg -framerate "$FPS" \
            -i "${IMG_DIR}/tiles_comps/${VAR1}-${VAR2}${VAR3:+-$VAR3}-%04d.png" \
            -hide_banner \
            -loglevel error \
            -c:v libx264 \
            -crf 23 \
            -pix_fmt yuv420p \
            -preset slow \
            -an -movflags +faststart \
            "$MOVIE_OUTPUT.mp4"
        echo "========================================================="
        echo "Movie created!"
    fi
fi

