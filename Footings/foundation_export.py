import re
from datetime import datetime

TEMPLATE_PATH = 'Footings/STK_PadFoundation_template.txt'
OUTPUT_PATH = ''

def create_mem_file(foundation_params):
    # Read the template file
    with open(TEMPLATE_PATH, 'r') as file:
        template_content = file.read()

    # Define multipliers for each parameter
    multipliers = {
        'L': 1000,
        'W': 1000,
        'D': 1000,
        'ColL': 1000,
        'ColW': 1000,
        'Pdl': 1,
        'Pll': 1,
        'BP': 1,
        'fc': 1,
        'Cvr': 1000,
        'barL': 1000,
        'CTSL': 1000,
        'barW': 1000,
        'CTSW': 1000
    }

    # Get current date and time
    current_datetime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Replace the placeholders with actual values
    for key, value in foundation_params.items():
        placeholder = f'[{key}]'
        if key in multipliers:
            # Apply multiplier and convert to integer
            value = int(value * multipliers[key])
        elif key == 'date':
            value = current_datetime
        template_content = template_content.replace(placeholder, str(value))

    # Replace any remaining placeholders with empty strings
    template_content = re.sub(r'\[.*?\]', '', template_content)

    # Write the modified content to the output file
    with open(OUTPUT_PATH+f'{foundation_params["foundation_ID"]}.mem', 'w') as file:
        file.write(template_content)
