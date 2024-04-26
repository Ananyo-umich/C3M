import yaml
import networkx as nx
import matplotlib.pyplot as plt
import re


def load_cantera_yaml(file_path):
    with open(file_path, "r") as file:
        return yaml.safe_load(file)


def create_reaction_graph(chemistry_data):
    G = nx.DiGraph()
    id = 0
    pattern = "[A-Za-z0-9_]+"
    for reaction in chemistry_data["reactions"]:
        print(reaction)
        reaction_id = id
        equation = reaction["equation"]
        reactants, products = map(str.strip, re.split("=>", equation))
        reactants = reactants.split(" + ")
        products = products.split(" + ")
        id = id + 1
        for reactant in reactants:
            reactant = re.findall(pattern, reactant)
            if reactant[-1] != "M":
                G.add_edge(reactant[-1], equation, type="reactant")

        for product in products:
            product = re.findall(pattern, product)
            if product[-1] != "M":
                G.add_edge(equation, product[-1], type="product")

    return G


def plot_reaction_graph(G):
    pos = nx.spring_layout(G)
    edge_labels = {(u, v): data["type"] for u, v, data in G.edges(data=True)}
    nx.draw(
        G,
        pos,
        with_labels=True,
        node_size=800,
        node_color="skyblue",
        font_size=8,
        font_weight="bold",
    )
    nx.draw_networkx_edge_labels(
        G, pos, edge_labels=edge_labels, font_size=8, font_color="red"
    )
    plt.title("Chemical Reaction Network")
    plt.show()


if __name__ == "__main__":
    cantera_yaml_file = "../tests/VenusBoxSOx.yaml"
    chemistry_data = load_cantera_yaml(cantera_yaml_file)
    reaction_graph = create_reaction_graph(chemistry_data)
    plot_reaction_graph(reaction_graph)
