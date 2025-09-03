import yaml


def upload_yaml(config_name):
    with open(f"./configs/{config_name}.yaml", "rb") as fp:
        config = yaml.load(fp, Loader=yaml.SafeLoader)
    return config
