import fitz

filepath = 'C:/Users/Adam/Desktop/temp/precastshopdrawings.pdf'

text = ''
with fitz.open(filepath) as doc:
    for page in doc:
        print("PAGE: "+str(page)+" *******************************")
        text = page.get_text()
        print(text)