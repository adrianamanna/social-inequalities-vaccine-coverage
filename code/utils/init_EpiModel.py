import yaml
import pickle


# .  INIT IFR
def init_IFR(levs, ifr_type, model_type):
    if ifr_type == "age":
        with open("./configs/IFR_age.yaml", "rb") as fp:
            IFR_age = yaml.load(fp, Loader=yaml.SafeLoader)
        if model_type == "Cij":
            IFR = IFR_age
        elif model_type == "Gab":
            IFR = {key: IFR_age[list(key)[0]] for key in levs}
    else:
        raise ValueError("Invalid `ifr_type`")
    return IFR


# .  UPLOADING DAILY VACCINATION DATA
def init_Omega(vax_scenario, nw, country):
    #  a. vaccination by age phases [as in Covid-19]
    if vax_scenario.startswith("age_structure"):
        with open(f"./data/vax_scenarios/Omega_t_age_{nw}.pkl", "rb") as f:
            Omega_t = pickle.load(f)
    else:
        raise ValueError("Invalid `vax_scenario`")
    return Omega_t
