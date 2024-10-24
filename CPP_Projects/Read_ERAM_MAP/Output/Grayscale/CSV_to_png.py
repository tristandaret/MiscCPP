import numpy as np
import pandas as pd
from PIL import Image

def csv_to_png(csv_file, output_image, image_format='PNG'):
    # Read the CSV file
    df = pd.read_csv(csv_file, header=None)
    
    # Convert the DataFrame to a numpy array
    data = df.to_numpy()
    
    # Convert to uint8
    data = data.astype(np.uint8)
    
    # Create an image from the numpy array
    image = Image.fromarray(data, mode='L')

    # Resize the image by the specified scale factor
    scale_factor = 25
    original_size = image.size
    new_size = (original_size[0] * scale_factor, original_size[1] * scale_factor)
    image = image.resize(new_size, Image.NEAREST)  # Use NEAREST to avoid introducing new colors
    
    # Save the image
    image.save(output_image, format=image_format)

# Do it for the 34 ERAMs
eram_indices = [1, 2, 3, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 23, 24, 26, 28, 29, 30, 36, 37, 38, 39, 41, 42, 43, 44, 45, 46, 47]

for i in range(len(eram_indices)):
    csv_file = f'Output/CSV/ERAM_{eram_indices[i]}_RC_grayscale.csv'
    output_image = f'Output/CSV/ERAM_{eram_indices[i]}_RC_grayscale_larger.png'
    csv_to_png(csv_file, output_image, 'PNG')
