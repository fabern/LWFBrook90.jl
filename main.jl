# fabian.bernhard@wsl.ch, 2021-01-02

using LWFBrook90Julia
greet()



input_defs = input_definitions("folder","prefix")
input_defs


TravisCI


# using PkgTemplates
# t = Template(;
#             user = "fabern",
#             license = "GPL-3.0+",
#             authors = ["Fabian Bernhard"],
            
#             plugins = [TravisCI(), Codecov(), BlueStyleBadge(),
#                         Documenter{TravisCI}()
#             ])

# generate("LWF", t)