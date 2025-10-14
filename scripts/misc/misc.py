import re
def to_snake_case(name):
    # Convert to lowercase, replace spaces and parentheses with underscores, remove multiple underscores
    name = name.lower()
    name = re.sub(r"[^\w]+", "_", name)
    name = re.sub(r"_+", "_", name)
    return name.strip("_")