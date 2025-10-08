import json
def json_dumps(data, sort=False, inline_limit=2, indent=0):
    """Serializes `data` into a JSON-formatted string.
     
    Serialization is performed with customized formatting:
    - Dictionaries are printed with each key/value pair on its own line.
    - Lists are printed inline if all their elements are simple types (str, int, float, bool, None).
      If any element is a list or a dict, the list is printed in a multiline format if such element
      is longer than inline_limit.
    - Indentation is done using tabs.
    
    Args:
        data (any): The data to serialize.
        sort (bool): If True, sort dictionary keys.
        inline_limit (int): Threshold on sublist elements size for printing inline.
        indent (int): Current nesting level.
    
    Returns:
        str: A JSON-formatted string.
    """
    sp = "\t" * indent  # current indentation using tabs

    if isinstance(data, dict):
        if sort:
            data = sorted(data.items())
        else:
            data = data.items()
        items = []
        for key, value in data:
            # Increase the indent level for nested content.
            item_str = f'{sp}\t{json.dumps(key)}: {json_dumps(value, indent=indent + 1, inline_limit=inline_limit)}'
            items.append(item_str)
        return '{\n' + ',\n'.join(items) + f'\n{sp}}}'
    
    elif isinstance(data, list):
        # Check if any element in the list is complex (a dict or a list)
        if any(isinstance(item, dict) or isinstance(item, list) and len(item) > inline_limit for item in data):
            items = []
            for item in data:
                # Increase the indent level for list items.
                item_str = json_dumps(item, indent=indent + 1, inline_limit=inline_limit)
                items.append("\t" * (indent + 1) + item_str)
            return "[\n" + ",\n".join(items) + "\n" + sp + "]"
        else:
            # Inline formatting for simple lists.
            items = [json_dumps(item, indent=0, inline_limit=inline_limit) for item in data]
            return '[' + ', '.join(items) + ']'
    
    else:
        # For primitive types, use the standard json.dumps for correct formatting.
        return json.dumps(data)