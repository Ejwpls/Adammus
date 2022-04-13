import os

pathname = '//hera-fs-01/Synergy\Projects/18/18129 Subiaco Pavilion/SHOP DRAWINGS/OUT/200929 - VERTICAL TRANSPORT - ROBERTS'

file_names = os.listdir(pathname)
suffix = "_HERA REVIEW"

new_names = file_names[:-len(".pdf")]

#[sub + suffix for sub in file_names]

print(new_names)




