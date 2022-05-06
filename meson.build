project('MinCostFlow','cpp',
    default_options: [
        'warning_level=3',
        'debug=true',
        'optimization=3',
        'c_std=c++17'])

subdir('include')

libmincostflow = library('mincostflow',
    my_headers,
    include_directories: my_includedir,
    install: true)

install_headers(my_headers,subdir: 'mincostflow')

#subdir('examples')
