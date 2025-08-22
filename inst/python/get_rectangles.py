from PIL import Image
import numpy as np
import cv2
#import csv
#import matplotlib.pyplot as plt  # For plotting the image
import re
import argparse
import os

from PIL import Image
import numpy as np
import cv2


def find_colored_contours_with_rectangles(image_path, target_color, color_name):
    """
    Identifies regions in a PNG image bounded by a specific color line and returns a series of rectangles
    as large as possible that fit within the contour without overlapping. Pixels that cannot fit into rectangles
    are returned individually.

    Parameters:
        image_path (str): Path to the input PNG image.
        target_color (tuple): RGB values to look for as the boundary color (e.g., (255, 0, 0) for red).
        color_name (str): Name of the color being detected.

    Returns:
        region_data (list): A list containing rectangles and individual pixel coordinates.
    """
    try:
        # Load the image using PIL
        image = Image.open(image_path).convert("RGB")
        image_array = np.array(image)  # Convert to a NumPy array for processing
        width, height = image.size
        
        # Create a mask for detecting the target color
        lower_bound = np.array([max(0,int(a) - 4) for a in target_color], dtype=np.uint8)
        upper_bound = np.array([min(int(a) + 5,255) for a in target_color], dtype=np.uint8)

        # Mask the image to extract pixels matching the target color
        color_mask = cv2.inRange(image_array, lower_bound, upper_bound)

        # Find contours from the masked image
        contours, _ = cv2.findContours(color_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

        # Prepare CSV output with region information
        region_data = []
        print('Found ' + str(len(contours)) + ' contours')

        for region_id, contour in enumerate(contours, start=1):  # Start region numbering from 1
            # Get all pixels within the contour
            mask = np.zeros_like(color_mask)
            cv2.drawContours(mask, [contour], -1, 255, thickness=cv2.FILLED)
            contained_pixels = np.array(np.where(mask == 255)).T  # Get pixel coordinates

            # Initialize grid representation of the contour
            contour_grid = np.zeros_like(mask)
            for pix in contained_pixels:
                contour_grid[pix[0], pix[1]] = 1

            # Extract rectangles from the contour grid
            rectangles = extract_rectangles(contour_grid, contour_grid.shape)

            # Add rectangles to region data
            for rect in rectangles:
                x_min, y_min, x_max, y_max = rect
                region_data.append({
                    "Color": color_name,
                    "Region_ID": region_id,
                    "x_pixels": y_min,
                    "y_pixels": height - x_max,
                    "Pixel_Count": (x_max - x_min) * (y_max - y_min),
                    "x_origin": 0,
                    "y_origin": 0,
                    "x_max": width,
                    "y_max": height,
                    "spot_height": x_max - x_min,
                    "spot_width": y_max - y_min,
                })

        return region_data

    except Exception as e:
        print(f"An error occurred: {e}")
        return None


def extract_rectangles(grid, grid_shape):
    """
    Extracts the largest possible rectangles from a binary grid representation of a contour.

    Args:
        grid (np.ndarray): Binary grid where "1" indicates pixels within the contour.
        grid_shape (tuple): Shape of the binary grid (rows, cols).

    Returns:
        list: A list of rectangle tuples (x_min, y_min, x_max, y_max).
    """
    rectangles = []
    visited = np.zeros(grid_shape, dtype=np.bool_)  # To keep track of visited pixels

    for i in range(grid_shape[0]):  # Iterate rows
        for j in range(grid_shape[1]):  # Iterate columns
            if grid[i, j] == 1 and not visited[i, j]:  # Found a starting pixel for a rectangle
                # Expand rectangle until we encounter a boundary in both dimensions
                x_min, y_min = i, j
                x_max, y_max = i, j

                # Expand downward
                while x_max + 1 < grid_shape[0] and np.all(grid[x_max + 1, y_min:y_max + 1] == 1) and np.all(~visited[x_max + 1, y_min:y_max + 1]):
                    x_max += 1

                # Expand rightward
                while y_max + 1 < grid_shape[1] and np.all(grid[x_min:x_max + 1, y_max + 1] == 1) and np.all(~visited[x_min:x_max + 1, y_max + 1]):
                    y_max += 1

                # Mark the pixels as visited
                visited[x_min:x_max + 1, y_min:y_max + 1] = True

                # Add the rectangle to the list
                rectangles.append((x_min, y_min, x_max + 1, y_max + 1))  # Rectangle bounds
    #print('Found '+str(len(rectangles))+' rectangles')
    return rectangles

def find_colored_contours(image_path, target_color, color_name):
    """
    Identifies regions in a PNG image bounded by a specific color line and saves the results.

    Parameters:
        image_path (str): Path to the input PNG image.
        target_color (tuple): RGB values to look for as the boundary color (e.g., (255, 0, 0) for red).
        output_plot_path (str): Path to save the plotted labeled image with detected contours.
        output_csv_path (str): Path to save the CSV file containing the pixel coordinates inside each region.

    """
    try:
        # Load the image using PIL
        image = Image.open(image_path).convert("RGB")
        image_array = np.array(image)  # Convert to a NumPy array for processing
        width, height = image.size
        # Create a mask for detecting the target color
        lower_bound = np.array([int(a)-2 for a in target_color], dtype=np.uint8)
        upper_bound = np.array([int(a)+2 for a in target_color], dtype=np.uint8)
        
        # Mask the image to extract pixels matching the target color
        color_mask = cv2.inRange(image_array, lower_bound, upper_bound)
        # Find contours from the masked image
        contours, _ = cv2.findContours(color_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        
        # Create a copy of the original image for drawing labeled contours
        labeled_image = image_array.copy()
        font_scale = 0.5
        font_thickness = 1
        
        # Prepare CSV output with region information
        region_data = []
        print('Found '+str(len(contours))+' contours')
        for region_id, contour in enumerate(contours, start=1):  # Start region numbering from 1
            # Label the contour on the image
            moments = cv2.moments(contour)
            if moments["m00"] != 0:  # Calculate centroid (avoid division by zero)
                center_x = int(moments["m10"] / moments["m00"])
                center_y = int(moments["m01"] / moments["m00"])
                cv2.putText(labeled_image, str(region_id), (center_x, center_y),
                            cv2.FONT_HERSHEY_SIMPLEX, font_scale, (0, 255, 0), font_thickness)
            
            # Get all pixels within the contour
            mask = np.zeros_like(color_mask)
            cv2.drawContours(mask, [contour], -1, 255, thickness=cv2.FILLED)
            contained_pixels = np.array(np.where(mask == 255)).T  # Get pixel coordinates

            # Add region data
            if len(contained_pixels) >2:
                ##now add in the pixels that didn't fit
                for pix in contained_pixels.tolist():
                  # if pix[0] > x+w | pix[0] < x | pix[1] >y | pix[1] < y-h:
                  region_data.append({
                        "Color": color_name,
                        "Region_ID": region_id,
                        "Pixel_Count": len(contained_pixels),
                        "x_pixels": pix[0],
                        "y_pixels": height-pix[1],
                        "x_origin": 0,
                        "y_origin": 0,
                        "x_max": width,
                        "y_max": height,
                        "spot_height": 1,
                        "spot_width": 1,
                  })
        return region_data

    except Exception as e:
        print(f"An error occurred: {e}")

def build_csv(image_path, output_csv_path, colorcodes):
    with open(output_csv_path,'w') as f:
        f.write("Color,Region_ID,Pixel_Count,x_pixels,y_pixels,x_origin,y_origin,x_max,y_max,spot_height,spot_width\n")
        for color in colorcodes:
            color_name, rgb = color.split(':')
            target_color = [int(a) for a in re.split(r',', rgb)]
            print(color_name, rgb)
            cdat = find_colored_contours_with_rectangles(image_path, target_color, color_name)
            for data in cdat:
                f.write(f"{data['Color']},{data['Region_ID']},{data['Pixel_Count']},{data['x_pixels']},{data['y_pixels']},{data['x_origin']},{data['y_origin']},{data['x_max']},{data['y_max']},{data['spot_height']},{data['spot_width']}\n")


        
if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('colorcodes',nargs = '+', help='List of colors to look for in the format of color:r,g,b')
    parser.add_argument('--image', help = 'Image file to process')

    args = parser.parse_args()
    image_path = args.image
    output_csv_path = re.sub('.png','_coords.csv',os.path.basename(image_path))

    build_csv(image_path, output_csv_path, args.colorcodes)
    

#for RA data:
#python process_image.py --image patient_031.png Blue:0,176,240 Yellow:255,255,0 Green:0,176,80

