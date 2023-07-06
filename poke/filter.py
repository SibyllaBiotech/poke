

class PocketFilter:
    def __init__(self, pocket):
        self.pocket = pocket

    def __hasitem__(self, key: str):
        return False

    def __getitem__(self, key: str):
        try:
            try:
                return getattr(self, key)
            except AttributeError:
                return getattr(self.pocket, key)
        except:
            raise KeyError(key)
    
    def evaluate(expr, pocket):
        return eval(expr, globals(), PocketFilter(pocket))
    

