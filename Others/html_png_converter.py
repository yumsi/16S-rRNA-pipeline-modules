import imgkit
import os
from htmlwebshot import WebShot, Config
from html2image import Html2Image

# biotaviz_sankey_prepfile-gender-Female.html
def main():
    html_to_png()


def html_to_png():
    shot = WebShot()
    # shot.flags = ["--quiet", "--enable-javascript", "--no-stop-slow-scripts", "--javascript-delay 1000", "--enable-local-file-access",
    #               "--run-script javascript: function(el) {d3.select(el).selectAll(".node text").filter(d => d.target).attr("x", -16).attr("text-anchor", "end");}"]
    shot.create_pic(html="biotaviz_sankey_prepfile-gender-Female.html", output="output.jpg")


main()
