from ...imports import *

__all__ = ["help"]


def help(self):
    """
    Print a quick reference of key actions available for this `Rainbow`.
    """
    print(
        textwrap.dedent(
            """
    Hooray for you! You asked for help on what you can do
    with this 🌈 object. Here's a quick reference of a few
    available options for things to try."""
        )
    )

    base_directory = pkg_resources.resource_filename("chromatic", "rainbows")
    descriptions_files = []
    for level in ["*", "*/*"]:
        descriptions_files += glob.glob(
            os.path.join(base_directory, level, "descriptions.txt")
        )
    categories = [
        d.replace(base_directory + "/", "").replace("/descriptions.txt", "")
        for d in descriptions_files
    ]
    for i in np.argsort(categories):
        c, d = categories[i], descriptions_files[i]
        header = (
            "\n" + "-" * (len(c) + 4) + "\n" + f"| {c} |\n" + "-" * (len(c) + 4) + "\n"
        )

        table = ascii.read(d)
        items = []
        for row in table:
            name = row["name"]
            if hasattr(self, name) or (name in ["+-*/", "[:,:]"]):
                if name in "+-*/":
                    function_call = f"{name}"
                else:
                    function_call = f".{name}()"

                item = (
                    f"{row['cartoon']} | {function_call:<28} \n   {row['description']}"
                )
                items.append(item)
        if len(items) > 0:
            print(header)
            print("\n".join(items))
