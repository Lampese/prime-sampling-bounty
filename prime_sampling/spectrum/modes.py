from importlib import import_module

MODES = {
    'single4': '.single4',
    'pair4': '.pair4',
    'six4': '.six4',
    'single6plus': '.single6plus',
    'single10': '.single10',
    'pair10loose': '.pair10loose',
    'pair10': '.pair10',
    'six10': '.six10',
}


def install_mode(mode, boxes):
    import_module(MODES[mode], __package__).install(boxes)
