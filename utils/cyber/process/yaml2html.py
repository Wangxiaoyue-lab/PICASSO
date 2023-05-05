import argparse
import yaml

# 解析命令行参数
parser = argparse.ArgumentParser()
parser.add_argument('input_file', help='Input YAML file')
parser.add_argument('output_file', help='Output HTML file')
#parser.add_argument('css_file', help='External CSS file')
args = parser.parse_args()

# 读取YAML文件
with open(args.input_file, 'r') as f:
    data = yaml.safe_load(f)

# 生成HTML文件
with open(args.output_file, 'w') as f:
    f.write('<html>\n')
    f.write('<head>\n')
    f.write('<style>\n')
    f.write('body { margin: 0; padding: 0; background-color: #D1EED1 }\n')
    f.write(
        '.container { background-color: white; width: 800px; margin: 2 auto; border-left: 1px solid #ccc; border-right: 1px solid #ccc; padding: 16px; }\n')
    f.write('h1, h2 { color: #333; }\n')
    f.write('h1 { font-size: 24px; margin: 16px 0; }\n')
    f.write('h2 { font-size: 18px; margin: 16px 0; }\n')
    f.write('p, li { color: #555;  font-size: 16px; line-height: 1.5; }\n')
    f.write('ul { margin: 0; padding: 0 }\n')
    f.write('li { list-style-type: none; }\n')
    f.write('</style>\n')
    f.write('<title>{}</title>\n'.format(data['Experiment_name']))
    #f.write('<link rel="stylesheet" type="text/css" href="{}">\n'.format(args.css_file))
    f.write('</head>\n')
    f.write('<body>\n')
    f.write('<div class="container">\n')
    f.write('<div style="background-color: green; text-align: center;height: 103px">\n')
    f.write(
        '<h1 style="font-family: \'Brush Script MT\', cursive; color: white; font-size: 78px">WangLab生信试验记录</h1>\n')
    f.write('</div>\n')
    f.write('<h1>实验名称:{}</h1>\n'.format(data['Experiment_name']))
    f.write('<p>实验单位：{}</p >\n'.format(data['Experiment_situation']['Experiment_location']))
    f.write('<p>实验者：{}</p >\n'.format(data['Experiment_situation']['Experiment_operater']))
    f.write('<p>实验日期：{}</p >\n'.format(data['Experiment_situation']['Experiment_date']))
    f.write('<p>实验平台：{}</p >\n'.format(data['Experiment_situation']['Experiment_platform']))
    f.write('<h1>实验目的</h1>\n')
    #f.write('<p>{}</p >\n'.format(data['实验目的']))
    f.write('<h1>实验方法与材料</h1>\n')
    f.write('<ul>\n')
    # for item in data['实验方法与材料']:
    #    f.write('<li>{}</li>\n'.format(item))
    f.write('</ul>\n')
    f.write('<h2>实验数据</h2>\n')
    f.write('<p>{}</p >\n'.format(data['experiment_methods_and_material']['Experiment_data']))
    f.write('<h2>实验软件</h2>\n')
    f.write('<ul>\n')
    for software in data['experiment_methods_and_material']['Experiment_software']:
        f.write('<li>{}</li>\n'.format(software))
    f.write('</ul>\n')
    f.write('<h2>补充</h2>\n')
    f.write('<p>{}</p >\n'.format(data['experiment_methods_and_material']['Supplement']))
    f.write('<h1>实验过程</h1>\n')
    f.write('<p>{}</p >\n'.format(data['Experiment_process']))
    f.write('<h1>实验结果</h1>\n')
    f.write('<p>{}</p >\n'.format(data['Experiment_result']))
    f.write('<h1>实验讨论</h1>\n')
    f.write('<p>{}</p >\n'.format(data['Experiment_discussion']))
    f.write('</body>\n')
    f.write('</div>\n')
