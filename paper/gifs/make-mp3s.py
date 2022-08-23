from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw 
from PIL import ImageOps
import os
import moviepy.editor as mp


FIG_DIRECTORY = "/Users/jtd/Documents/research/12.420 project/code/density/figs"
IMAGE_LAYOUTS = {
    "figure-5-animated": [
        ("Asymmetric reference asteroid",),
        ("Finite element model", "avg-asym"),
        ("Lumpy model", "asym-lump"),
        ("Symmetric reference asteroid",),
        ("Finite element model", "avg-sym"),
        ("Lumpy model", "sym-lump"),
    ],
    "figure-6-animated": [
        ("Finite element model", "avg-move-1.5"),
        ("Lumpy model", "move-1.5-lump"),
    ],
    "figure-7-animated": [
        ("Finite element model", "avg-sph-3"),
        ("Lumpy model", "sph-3-lump"),
    ],
    "figure-8-animated": [
        ("Finite element model", "avg-double"),
        ("Lumpy model", "double-lump"),
    ],
}
RESIZE_FRAC = 0.5
HORIZONTAL_BUFFER = 24
VERTICAL_TEXT_WIDTH = 24
font = ImageFont.truetype("cmunrm.ttf", 24)

def make_gifs():
    for save_name, layout in IMAGE_LAYOUTS.items():
        print(save_name)
        images = []
        first_horizontal_text = True
        offset_y = 0
        for row in layout:
            if len(row) == 1:
                if not first_horizontal_text:
                    offset_y += HORIZONTAL_BUFFER # Buffer before text
                first_horizontal_text = False
                offset_y += HORIZONTAL_BUFFER # Offset for horizontal text
            if len(row) == 2:
                images.append([
                    Image.open(f'{FIG_DIRECTORY}/{row[1]}/fe-d.gif'),
                    Image.open(f'{FIG_DIRECTORY}/{row[1]}/fe-s.gif'),
                    Image.open(f'{FIG_DIRECTORY}/{row[1]}/fe-u.gif'),
                    Image.open(f'{FIG_DIRECTORY}/{row[1]}/fe-r.gif')
                ])
        full_offset_y = offset_y
        
        num_frames = None
        single_size = None
        duration = None
        for line in images:
            for im in line:
                if num_frames is None:
                    num_frames = im.n_frames
                    single_size = im.size
                    duration = im.info['duration']
                if im.n_frames != num_frames:
                    num_frames = min(num_frames, im.n_frames)
                if im.size != single_size:
                    raise Exception(f"Got image with {im.size} size, but expected {single_size}")
                if im.info['duration'] != duration:
                    duration = max(duration, im.info['duration'])

        # Make the images
        single_size = (int(single_size[0] * RESIZE_FRAC), int(single_size[1] * RESIZE_FRAC))
        frame_size = [
            single_size[0] * 4 + VERTICAL_TEXT_WIDTH,
            single_size[1] * len(images) + full_offset_y]
        frames = []
        for frame_index in range(num_frames):
            frame = Image.new('RGB', frame_size, color=(255, 255, 255))
            draw = ImageDraw.Draw(frame)

            # Images first
            image_row = 0
            offset_y = 0
            first_horizontal_text = True
            for row in layout:
                if len(row) == 1:
                    # Horizontal text
                    if not first_horizontal_text:
                        offset_y += HORIZONTAL_BUFFER # Buffer before text
                    first_horizontal_text = False
                    offset_y += HORIZONTAL_BUFFER
                else:
                    gif_frame = int(frame_index / num_frames * images[image_row][0].n_frames)
                    for image_index in range(4):
                        images[image_row][image_index].seek(gif_frame)
                        frame.paste(images[image_row][image_index].resize(single_size),
                            (image_index * single_size[0] + VERTICAL_TEXT_WIDTH, image_row * single_size[1] + offset_y))
                    image_row += 1

            # Now the fonts
            image_row = 0
            offset_y = 0
            first_horizontal_text = True
            for row in layout:
                if len(row) == 1:
                    # Horizontal text
                    if not first_horizontal_text:
                        offset_y += HORIZONTAL_BUFFER # Buffer before text
                    first_horizontal_text = False

                    line_width = int(draw.textlength(row[0], font=font))
                    draw.text((
                        (frame_size[0] - line_width) / 2,
                        image_row * single_size[1] + offset_y), row[0], (0,0,0), font=font)

                    offset_y += HORIZONTAL_BUFFER
                else:
                    # Write text
                    line_width = int(draw.textlength(row[0], font=font))
                    text = Image.new('RGBA', (line_width,48))
                    text_draw = ImageDraw.Draw(text)
                    text_draw.text((0, 0), row[0], (0,0,0), font=font)
                    rotated=text.rotate(90,  expand=1)
                    frame.paste(rotated,(
                        0,
                        image_row * single_size[1] + offset_y + (single_size[1] - line_width) // 2), rotated)
                    image_row += 1

            frames.append(frame)
        
        # Save as gif
        print("Saving gif")
        frames[0].save(fp=f"{save_name}.gif", format='GIF', append_images=frames,
                save_all=True, duration=duration, loop=0)

        # gif to mp4
        print("Converting to mp4")
        clip = mp.VideoFileClip(f"{save_name}.gif")
        clip.write_videofile(f"{save_name}.mp4")
        # os.system(f"ffmpeg -f gif -i \"{save_name}.gif\" \"{save_name}.mp4\"")

def make_trifecta():
    images = [
        Image.open(f'{FIG_DIRECTORY}/trifecta/fe-d.gif'),
        Image.open(f'{FIG_DIRECTORY}/trifecta/fe-t.gif'),
        Image.open(f'{FIG_DIRECTORY}/trifecta/fe-u.gif'),
    ]
    
    num_frames = None
    single_size = None
    duration = None
    for im in images:
        if num_frames is None:
            num_frames = im.n_frames
            single_size = im.size
            duration = im.info['duration']
        if im.n_frames != num_frames:
            num_frames = min(num_frames, im.n_frames)
        if im.size != single_size:
            raise Exception(f"Got image with {im.size} size, but expected {single_size}")
        if im.info['duration'] != duration:
            duration = max(duration, im.info['duration'])

    # Make the images
    single_size = (int(single_size[0] * RESIZE_FRAC), int(single_size[1] * RESIZE_FRAC))
    frame_size = [
        single_size[0] * 2,
        single_size[1] * 2]
    frames = []
    for frame_index in range(num_frames):
        frame = Image.new('RGB', frame_size, color=(255, 255, 255))

        # Images first
        for i in range(len(images)):
            gif_frame = int(frame_index / num_frames * images[i].n_frames)
            images[i].seek(gif_frame)
            if i == 0:
                frame.paste(images[i].resize(single_size),
                    (0, 0))
            if i == 1:
                frame.paste(images[i].resize(single_size),
                    (single_size[0], 0))
            if i == 2:
                frame.paste(images[i].resize(single_size),
                    (single_size[0] // 2, single_size[1]))
        frames.append(frame)
    
    # Save as gif
    print("Saving gif")
    frames[0].save(fp=f"figure-4-animated.gif", format='GIF', append_images=frames,
            save_all=True, duration=duration, loop=0)

    # gif to mp4
    print("Converting to mp4")
    clip = mp.VideoFileClip(f"figure-4-animated.gif")
    clip.write_videofile(f"figure-4-animated.mp4")

if __name__ == "__main__":
    make_gifs()
    make_trifecta()