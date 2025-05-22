import argparse
import torch
from recbole.quick_start import run_recbole

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--model', '-m', type=str, default='FRF', help='模型名称')
    parser.add_argument('--dataset', '-d', type=str, default='ml-1m', help='数据集名称')
    parser.add_argument('--config_files', type=str, default='FreqDualDynamic.yaml', help='配置文件路径')
    parser.add_argument('--params', type=str, default=None, help='模型超参数')

    args, _ = parser.parse_known_args()

    # 配置参数
    config_file_list = args.config_files.strip().split(' ') if args.config_files else None
    
    # 命令行参数
    model_params = {}
    if args.params:
        param_list = args.params.strip().split(' ')
        for param in param_list:
            key, value = param.split('=')
            model_params[key] = value
            
    # 运行模型
    print(f'Running {args.model} on {args.dataset}')
    run_recbole(
        model=args.model,
        dataset=args.dataset,
        config_file_list=config_file_list,
        config_dict=model_params
    )