import os
import re
import textwrap
from jinja2 import Template
import yaml
import subprocess

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
WRAPPER_DIR = os.path.dirname(BASE_DIR)
OUTPUT_DIR = os.path.join(BASE_DIR, "wrappers")
META_OUTPUT_DIR = os.path.join(BASE_DIR, "meta-wrappers")
SCRIPTS = {"wrapper.py", "wrapper.R"}
BLACKLIST = {
    ".snakemake",
    ".circleci",
    "wercker.yml",
    "docs",
    "__pycache__",
    "environment.yaml",
    ".git",
    ".github",
    ".gitignore",
    "README.md",
    ".cache",
    "__init__.py",
    "test_wrappers.py",
    ".pytest_cache",
    ".gitkeep",
} | SCRIPTS

try:
    TAG = subprocess.check_output(["git", "describe", "--tags"]).decode().strip()
except subprocess.CalledProcessError:
    TAG = "No_Tags"

DISCIPLINES = [
    os.path.basename(x) for x in os.listdir(WRAPPER_DIR) if all(
        [x not in BLACKLIST, os.path.isdir(os.path.join(WRAPPER_DIR, x)), x!="meta"])
    ]
META_DISCIPLINES = [
    os.path.basename(x) for x in os.listdir(
        os.path.join(WRAPPER_DIR, "meta")) if all(
            [x not in BLACKLIST, os.path.isdir(os.path.join(WRAPPER_DIR, "meta", x))]
        )]
disciplines_mask = {
    "bio": {"name": "Bio", "description": "bioinformatics"},
    "geo": {"name": "Geo", "description": "geospatial"},
    "phys": {"name": "Phys", "description": "physics"},
    "utils": {"name": "Utils", "description": "utility"}
}

with open(os.path.join(BASE_DIR, "_templates", "tool.rst")) as f:
    TOOL_TEMPLATE = Template(f.read())

with open(os.path.join(BASE_DIR, "_templates", "wrapper.rst")) as f:
    TEMPLATE_WRAPPER = Template(f.read(), trim_blocks=True, lstrip_blocks=True)

with open(os.path.join(BASE_DIR, "_templates", "meta_wrapper.rst")) as f:
    TEMPLATE_META = Template(f.read())

with open(os.path.join(BASE_DIR, "_templates", "discipline.rst")) as f:
    TEMPLATE_DISCIPLINE = Template(f.read())

with open(os.path.join(BASE_DIR, "_templates", "meta_discipline.rst")) as f:
    TEMPLATE_META_DISCIPLINE = Template(f.read())

def get_tool_dir(tool, discipline=None):
    outdir = os.path.join(OUTPUT_DIR, discipline, tool)
    os.makedirs(outdir, exist_ok=True)
    return outdir


def render_tool(tool, discipline, subcmds):
    if "" in subcmds:
        raise NotImplementedError(
            'You are trying to render the tool "',
            tool,
            '" for subcommands "',
            subcmds,
            "\". The empty subcommand '' indicates that both the main command, as well ",
            "as the subcommand(s) have wrappers. This case requires a Template hybrid ",
            "between wrapper.rst and tool.rst that has not been implemented.",
        )
    with open(os.path.join(get_tool_dir(tool, discipline) + ".rst"), "w") as f:
        f.write(TOOL_TEMPLATE.render(name=tool))


def render_snakefile(path):
    if not os.path.isfile(os.path.join(path, "test", "Snakefile")):
        return "test/Snakefile not found."
    
    with open(os.path.join(path, "test", "Snakefile")) as snakefile:
        lines = filter(lambda line: re.search(r"# ?\[hide\]", line) is None, snakefile)
        snakefile = textwrap.indent(
            "\n".join(l.rstrip() for l in lines).strip(), "    "
        ).replace("master", TAG)
        return snakefile


def render_wrapper(path, target, wrapper_id):
    print("rendering", path)
    with open(os.path.join(path, "meta.yaml")) as meta:
        meta = yaml.load(meta, Loader=yaml.BaseLoader)
    if "blacklisted" not in meta:
        meta["blacklisted"] = None

    envpath = os.path.join(path, "environment.yaml")
    if os.path.exists(envpath):
        with open(envpath) as env:
            env = yaml.load(env, Loader=yaml.BaseLoader)
            if env:
                # pkgs = env["dependencies"]
                pkgs = env.get("dependencies", [])
            else:
                pkgs = []
    else:
        pkgs = []

    snakefile = render_snakefile(path)

    wrapper = os.path.join(path, "wrapper.py")
    wrapper_lang = "python"
    if not os.path.exists(wrapper):
        wrapper = os.path.join(path, "wrapper.R")
        wrapper_lang = "R"
    with open(wrapper) as wrapper:
        wrapper = textwrap.indent(wrapper.read(), "    ")

    name = meta["name"].replace(" ", "_") + ".rst"
    os.makedirs(os.path.dirname(target), exist_ok=True)
    with open(target, "w") as readme:
        rst = TEMPLATE_WRAPPER.render(
            snakefile=snakefile,
            wrapper=wrapper,
            wrapper_lang=wrapper_lang,
            pkgs=pkgs,
            id=wrapper_id,
            wrapper_path=os.path.relpath(path, WRAPPER_DIR),
            **meta,
        )
        readme.write(rst)
    return name


def render_meta(path, target):
    print("rendering", path)
    with open(os.path.join(path, "meta.yaml")) as meta:
        meta = yaml.load(meta, Loader=yaml.BaseLoader)
    wrapperpath = os.path.join(path, "used_wrappers.yaml")
    if os.path.exists(wrapperpath):
        with open(wrapperpath) as env:
            env = yaml.load(env, Loader=yaml.BaseLoader)
            used_wrappers = env["wrappers"]
    else:
        used_wrappers = []
    snakefile = render_snakefile(path)

    name = meta["name"].replace(" ", "_") + ".rst"
    os.makedirs(os.path.dirname(target), exist_ok=True)
    with open(target, "w") as readme:
        rst = TEMPLATE_META.render(
            snakefile=snakefile,
            usedwrappers=used_wrappers,
            **meta,
        )
        readme.write(rst)
    return name


def setup(*args):
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    for discipline in DISCIPLINES:
        # generate discipline index
        with open(os.path.join(OUTPUT_DIR, discipline + ".rst"), "w") as f:
            f.write(TEMPLATE_DISCIPLINE.render(
                name=disciplines_mask[discipline]['name'],
                description=disciplines_mask[discipline]['description']))
        os.makedirs(os.path.join(OUTPUT_DIR, discipline), exist_ok=True)
        for tool in os.listdir(os.path.join(WRAPPER_DIR, discipline)):
            if tool in BLACKLIST:
                continue
            subcmds = []
            max_depth = 1
            path = os.path.join(WRAPPER_DIR, discipline, tool)
            for dirpath, dirnames, filenames in os.walk(path):
                if any(f in SCRIPTS for f in filenames):
                    # split up paths including SCRIPTS on the displine and strip tool from start
                    tool_levels = dirpath.rpartition(os.sep.join(["", discipline, ""]))[
                        2
                    ].split(sep=os.sep)
                    wrapper_id = os.path.join(discipline, *tool_levels)
                    tool_levels.remove(tool)
                    if len(tool_levels) == 0:
                        continue
                    if any(l in BLACKLIST for l in tool_levels):
                        continue
                    subcommand = os.sep.join(tool_levels)
                    subcmds.append(subcommand)
                    # watermark max_depth for current tool
                    max_depth = max(max_depth, len(tool_levels))
                    render_wrapper(
                        dirpath,
                        os.path.join(get_tool_dir(tool, discipline), subcommand + ".rst"),
                        wrapper_id,
                    )
            if len(subcmds) == 0:
                # no subcommands, render a wrapper for the main tool command
                render_wrapper(
                    path,
                    os.path.join(OUTPUT_DIR, discipline, tool + ".rst"),
                    os.path.join(discipline, tool),
                )
            else:
                # subcommands found (and rendered above), render the tool's table of content
                render_tool(tool, discipline, subcmds)

    os.makedirs(META_OUTPUT_DIR, exist_ok=True)
    for discipline in META_DISCIPLINES:
        # generate meta-discipline index
        with open(os.path.join(META_OUTPUT_DIR, discipline + ".rst"), "w") as f:
            f.write(TEMPLATE_META_DISCIPLINE.render(**disciplines_mask[discipline]))
        for subwf in os.listdir(os.path.join(WRAPPER_DIR, "meta", discipline)):
            if subwf in BLACKLIST:
                continue
            path = os.path.join(WRAPPER_DIR, "meta", discipline, subwf)
            render_meta(path, os.path.join(META_OUTPUT_DIR, discipline, subwf + ".rst"))


if __name__ == "__main__":
    setup()
