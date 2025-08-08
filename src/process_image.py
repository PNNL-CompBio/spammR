from PIL import Image
import numpy as np
import cv2
import csv
import matplotlib.pyplot as plt  # For plotting the image
import re
import argparse
import os

from PIL import Image
import cv2
import numpy as np
import matplotlib.pyplot as plt

def find_colored_contours(image_path, target_color, output_plot_path, output_csv_path):
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
        
        # Create a mask for detecting the target color
        lower_bound = np.array([a-2 for a in target_color], dtype=np.uint8)
        upper_bound = np.array([a+2 for a in target_color], dtype=np.uint8)
        
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
            region_data.append({
                "Region_ID": region_id,
                "Pixel_Count": len(contained_pixels),
                "Pixels": contained_pixels.tolist()
            })
        
        # Save region data to a CSV file
        with open(output_csv_path, 'w') as f:
            f.write("Region_ID,Pixel_Count,Pixels\n")
            for data in region_data:
                f.write(f"{data['Region_ID']},{data['Pixel_Count']},{data['Pixels']}\n")
        
        # Plot the labeled image with regions
        plt.figure(figsize=(12, 12))
        plt.imshow(cv2.cvtColor(labeled_image, cv2.COLOR_BGR2RGB))
        plt.title("Labeled Contours")
        plt.axis("off")
        plt.savefig(output_plot_path)
        plt.show()
        print("Plotted labeled image with contours saved.")
        print("Output CSV saved with region data.")
    except Exception as e:
        print(f"An error occurred: {e}")

# Example usage:
#image_path = "test_image.png"
#target_color = (255, 0, 0)  # Example target color (red)
#output_plot_path = "labeled_contours.png"
#output_csv_path = "contour_regions.csv"



if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--image', help = 'Image file to process')
    parser.add_argument('--rgb', help = 'comma-delimited RGB color to track', default='7,173,236')
    parser.add_argument('--color',help='Name of color provided',default='blue')

    args = parser.parse_args()
    image_path = args.image
    target_color = [int(a) for a in re.split(r',', args.rgb)]
    color_name = args.color
    output_csv_path = re.sub('.png',target_color+'_coords.csv',os.path.basename(image_path))
    output_plot_path = re.sub('.png',target_color+'_output.png',os.path.basename(image_path))
    find_colored_contours(image_path, target_color, output_plot_path, output_csv_path)
#save_pixels_and_plot_contours(image_path, output_file,int(threshold))
# Example Usage:
#image_path = "example.png"  # Replace with the path to your PNG image

                         #output_file = "pixels_by_contour_with_region_id.csv"  # Path to the CSV file for output

