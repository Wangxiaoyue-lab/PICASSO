import argparse
import yaml

# 解析命令行参数
parser = argparse.ArgumentParser()
parser.add_argument('input_file', help='Input YAML file')
parser.add_argument('output_file', help='Output HTML file')
parser.add_argument('css_file', help='External CSS file')
args = parser.parse_args()

# 读取YAML文件
with open(args.input_file, 'r') as f:
    data = yaml.safe_load(f)

# 生成HTML文件
with open(args.output_file, 'w') as f:
    f.write('<html>\n')
    f.write('<head>\n')
    f.write('<title>{}</title>\n'.format(data['实验名称']))
    f.write('<link rel="stylesheet" type="text/css" href="{}">\n'.format(args.css_file))
    f.write('</head>\n')
    f.write('<body>\n')
    f.write('<h1>{}</h1>\n'.format(data['实验名称']))
    f.write('<p>实验单位：{}</p >\n'.format(data['实验单位']))
    f.write('<p>实验者：{}</p >\n'.format(data['实验者']))
    f.write('<p>实验日期：{}</p >\n'.format(data['实验日期']))
    f.write('<p>实验平台：{}</p >\n'.format(data['实验平台']))
    f.write('<h2>实验目的</h2>\n')
    f.write('<p>{}</p >\n'.format(data['实验目的']))
    f.write('<h2>实验方法与材料</h2>\n')
    f.write('<ul>\n')
    for item in data['实验方法与材料']:
        f.write('<li>{}</li>\n'.format(item))
    f.write('</ul>\n')
    f.write('<h2>实验数据</h2>\n')
    f.write('<p>{}</p >\n'.format(data['实验数据']))
    f.write('<h2>实验软件</h2>\n')
    f.write('<ul>\n')
    for software in data['实验软件']:
        f.write('<li>{}</li>\n'.format(software))
    f.write('</ul>\n')
    f.write('<h2>补充</h2>\n')
    f.write('<p>{}</p >\n'.format(data['补充']))
    f.write('<h2>实验过程</h2>\n')
    f.write('<p>{}</p >\n'.format(data['实验过程']))
    f.write('<h2>实验结果</h2>\n')
    f.write('<p>{}</p >\n'.format(data['实验结果']))
    f.write('<h2>实验讨论</h2>\n')
    f.write('<p>{}</p >\n'.format(data['实验讨论']))
    f.write('</body>\n')