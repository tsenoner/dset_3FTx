import os
from PIL import Image, ImageDraw, ImageFont
import pandas as pd
import distinctipy
from pathlib import Path
import matplotlib.font_manager
import argparse

def generate_colors(num_colors, colorblind_type="Normal", seed=42):
    colors = distinctipy.get_colors(num_colors, colorblind_type=colorblind_type, rng=seed)
    colors = [[int(c * 255) for c in color] for color in colors]
    return colors

def create_legend(image_path, output_path, group_colors, font_size=12, shape_size=50, padding=10):
    im = Image.open(image_path)

    arial_font_path = matplotlib.font_manager.findfont("Arial")
    font = ImageFont.truetype(arial_font_path, font_size)
    draw = ImageDraw.Draw(im)

    legend_width = max([draw.textlength(g, font=font) for g in group_colors.keys()]) + shape_size + padding * 3
    legend_height = len(group_colors) * (shape_size + padding) + padding

    legend = Image.new("RGBA", (int(legend_width), int(legend_height)), (255, 255, 255, 0))
    legend_draw = ImageDraw.Draw(legend)

    for idx, (group, color) in enumerate(sorted(group_colors.items())):
        y_pos = idx * (shape_size + padding) + padding
        legend_draw.ellipse([padding, y_pos, padding + shape_size, y_pos + shape_size], fill=tuple(color))
        legend_draw.text((padding * 2 + shape_size, y_pos), group, font=font, fill=(0, 0, 0))

    legend_draw.rectangle([0, 0, int(legend_width) - 1, int(legend_height) - 1], outline=(0, 0, 0), width=1)
    im.paste(legend, (im.width - int(legend_width) - padding, padding), legend)
    im.save(output_path)

def add_legend_to_image(input_image, group_colors_csv=None, group_column=None, group_colors=None, output_image=None, font_size=12, shape_size=50):
    print("test")
    if group_colors is None:
        df = pd.read_csv(group_colors_csv)
        num_colors = len(df[group_column].unique())
        colors = generate_colors(num_colors)
        group_colors = {g: c for g, c in zip(df[group_column].unique(), colors)}
    if not output_image:
        output_image = input_image.with_stem(input_image.stem + "_legend")

    create_legend(input_image, output_image, group_colors, font_size=font_size, shape_size=shape_size)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add a legend to a CLANS plot")
    parser.add_argument("input_image", type=Path, help="Path to the input image")
    parser.add_argument("-o", "--output_image", type=Path, help="Path to the output image")
    parser.add_argument("-g", "--group_colors_csv", type=Path, help="Path to the CSV file with group colors")
    parser.add_argument("-c", "--group_column", help="Group column name in the CSV file")
    parser.add_argument("-f", "--font_size", type=int, default=12, help="Font size of the legend values")
    parser.add_argument("-s", "--shape_size", type=int, default=50, help="Shape size of the legend")
    args = parser.parse_args()
    add_legend_to_image(args.input_image, args.group_colors_csv, args.group_column, output_image=args.output_image, font_size=args.font_size, shape_size=args.shape_size)

