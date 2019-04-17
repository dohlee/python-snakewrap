def alwayslist(x):
    if not isinstance(x, list):
        return [x]
    return x