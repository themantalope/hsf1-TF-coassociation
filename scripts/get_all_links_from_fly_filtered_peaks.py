import re

link_pattern = re.compile('"((http|ftp)s?://.*?)"')
# url = "http://encode-ftp.s3.amazonaws.com/modENCODE_VS_ENCODE/Regulation/Fly/peakCalls/finalPk/"

# page = requests.get(url)

with open("finalPk.html", "r") as f:
    text = f.read()


all_links = re.findall(link_pattern, text)

good_links = [l[0] for l in all_links if l[1] == "http" and "narrowPeak.gz" in l[0]]
# print good_links

with open("fly_narrowpeaks.txt", "w") as f:
    f.write("\n".join(good_links))

