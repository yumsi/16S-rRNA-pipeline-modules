#!/usr/bin/python3

# Author: Harm Laurense
# Last edited: 01-11-2021
# Function: This script is used for creating .png files by making creenshots of sankey html files. These files are
# also adjusted by adding an overlay/text to show taxonomic ranks + title. These files are then merged into a grid (png)
# TODO: This script potentially needs altering for it to be implemented into nextflow.
# TODO: Add docstrings to every function
# TODO: Bugtesting / limittesting; how many files is too big for the grid to handle (for example)

import glob
from PIL import Image, ImageDraw, ImageFont  # PIL package itself is outdated. Now found inside the pillow package
from math import ceil
from selenium import webdriver


def main():
    screenshot_sankey()
    make_grid()


def draw_overlay(path_to_file, filename_notitle):
    img = Image.open(path_to_file + filename_notitle)
    font = ImageFont.truetype('C:\\Windows\\Fonts\\ARIALNB.ttf', 90)
    draw = ImageDraw.Draw(img)
    # Place filename as title for each sankey diagram
    title = filename_notitle.split(".png")[0]
    draw.text((100, 50), title, (0, 0, 0), font=font)

    # Place taxonomic rank overlay on each sankey diagram
    font = ImageFont.truetype('C:\\Windows\\Fonts\\ARIALNB.ttf', 50)
    overlay_domain = "Domain"
    draw.text((0, 3460), overlay_domain, (0, 0, 0), font=font)
    overlay_taxonomy = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]
    horizontal_spacing = 280
    correction = 0
    for taxonomic_rank in overlay_taxonomy:
        draw.text((horizontal_spacing + correction, 3465), taxonomic_rank, (0, 0, 0), font=font)
        horizontal_spacing += 835
        correction += 15
    img.save(path_to_file + filename_notitle)


def make_grid():
    frame_width = 4000
    images_per_row = 2
    padding = 10
    images = []
    filename_grid = "sankeyGrid"

    for filename in glob.glob("C:\\Users\\Harm_\\PycharmProjects\\GraduationTesting\\results - png sankey\\*.png"):
        images.append(filename)

    img_width, img_height = Image.open(images[0]).size
    sf = (frame_width - (images_per_row - 1) * padding) / (images_per_row * img_width)  # scaling factor
    scaled_img_width = ceil(img_width * sf)
    scaled_img_height = ceil(img_height * sf)
    number_of_rows = ceil(len(images) / images_per_row)
    frame_height = ceil(sf * img_height * number_of_rows)

    new_im = Image.new('RGB', (frame_width, frame_height))

    i, j = 0, 0
    for num, im in enumerate(images):
        if num % images_per_row == 0:
            i = 0
        im = Image.open(im)
        # Here I resize my opened image, so it is no bigger than 100,100
        im.thumbnail((scaled_img_width, scaled_img_height))
        # Iterate through a 4 by 4 grid with 100 spacing, to place my image
        y_cord = (j // images_per_row) * scaled_img_height
        new_im.paste(im, (i, y_cord))
        # print(i, y_cord)
        i = (i + scaled_img_width) + padding
        j += 1
    new_im.save(r"C:\Users\Harm_\PycharmProjects\GraduationTesting\\" + filename_grid + ".png",
                "PNG", quality=100, optimize=True, progressive=True)


def screenshot_sankey():
    files = []
    for filename in glob.glob("C:\\Users\\Harm_\\PycharmProjects\\GraduationTesting\\html files\\*.html"):
        files.append(filename)
    options = webdriver.ChromeOptions()
    options.headless = True
    driver = webdriver.Chrome(options=options)
    for html_file in files:
        html_file = html_file.split("\\")[-1]
        driver.get("file://C:/Users/Harm_/PycharmProjects/GraduationTesting/html files/" + html_file)
        S = lambda X: driver.execute_script('return document.body.parentNode.scroll' + X)
        driver.set_window_size(S('Width'), S('Height'))  # May need manual adjustment
        size = driver.get_window_size()
        print(size)
        html_file_name = html_file.split(".html")[0] + ".png"
        path = "C:\\Users\\Harm_\\PycharmProjects\\GraduationTesting\\results - png sankey\\"
        driver.find_element_by_tag_name('body').screenshot(path + html_file_name)
        draw_overlay(path, html_file_name)
    driver.quit()


main()
