import re

link_pattern = re.compile('"((http|ftp)s?://.*?)"')

with open("human_finalPk.html", "r") as f:
    text = f.read()

# print text

all_links = re.findall(link_pattern, text)
# print all_links

good_links = [l[0] for l in all_links if l[1] == "http" and ".gz" in l[0]]
# print good_links

with open("human_narrowpeaks.txt", "w") as f:
    f.write("\n".join(good_links))