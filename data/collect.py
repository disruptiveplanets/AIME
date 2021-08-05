import os, sys
from PIL import Image, ImageFont, ImageDraw

PATH = 'converge'
RUNS = 0
CORNER = 1
PARAMS = 2
REDCHI = 3
font = ImageFont.truetype("/usr/share/fonts/truetype/cmu/cmunrm.ttf", 20)

dirs = next(os.walk(PATH))[1]
for d in dirs:
    try:
        images = [Image.open(PATH+'/'+d+'/'+d+x) for x in ['-compare.png', '-corner.png', '-params.png', '-redchi.png']]
    except:
        continue
    widths, heights = zip(*(i.size for i in images))

    bi_height = max(heights[PARAMS], heights[REDCHI])
    bi_width = widths[PARAMS] + widths[REDCHI]

    im_width = widths[CORNER] + max(widths[RUNS], bi_width)
    im_height = max(heights[CORNER], heights[RUNS] + bi_height)

    new_im = Image.new('RGB', (im_width, im_height), color=(255, 255, 255))

    new_im.paste(images[CORNER], (0, int((im_height - heights[CORNER]) / 2)))
    new_im.paste(images[RUNS],
        (widths[CORNER] + int((im_width - widths[CORNER] - widths[RUNS]) / 2),
        int((im_height - heights[RUNS] - bi_height) / 2)))
    new_im.paste(images[PARAMS],
        (widths[CORNER] + int((im_width - widths[CORNER] - bi_width) / 2),
        heights[RUNS] + int((im_height - heights[RUNS] - bi_height) / 2)))
    new_im.paste(images[REDCHI],
        (widths[PARAMS] + widths[CORNER] + int((im_width - widths[CORNER] - bi_width) / 2),
        heights[RUNS] + int((im_height - heights[RUNS] - bi_height) / 2)))

    f = open(PATH+"/"+d+"/"+d+".dat", 'r')
    cadence = int(f.readline())
    impact_parameter = int(f.readline())
    speed = float(f.readline())
    spin = [float(x) for x in f.readline().split(',')]
    jlms = [float(x) for x in f.readline().split(',')]
    theta_true = [float(x) for x in f.readline().split(',')]
    theta_start = [float(x) for x in f.readline().split(',')]
    theta_spread = [float(x) for x in f.readline().split(',')]
    theta_high = [float(x) for x in f.readline().split(',')]
    theta_low = [float(x) for x in f.readline().split(',')]
    sigma = float(f.readline())
    f.close()

    draw = ImageDraw.Draw(new_im)
    text = """Cadence: {}
Speed: {}\t Impact parameter: {} Earth radii
Spin: {}
Jlms: {}
True theta: {}
Theta start: {}
Theta spread: {}
Theta high: {}
Theta low: {}
""".format(cadence, speed, impact_parameter, spin, jlms, theta_true, theta_start,
    theta_spread, theta_high, theta_low)

    draw.text((widths[CORNER]/2+40, 40), text, font=font, align="left", fill='black')
    new_im.save(PATH+'/'+d+'/'+d+"-all.png")

    for x in ['-compare.png', '-corner.png', '-params.png', '-redchi.png']:
        os.remove(PATH+'/'+d+'/'+d+x) 
