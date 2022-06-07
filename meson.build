project('MinCostFlow','cpp',
    default_options: [
        'warning_level=3',
        'debug=true',
        'optimization=3',
        'cpp_std=c++17'])

pkg=import('pkgconfig')
python3=dependency('python3')
boost=dependency('boost',modules: ['python3.7'])

subdir('include')
# subdir('src')

libmincostflow = library('mincostflow',
    my_headers,
    include_directories: my_includedir,
    install: true)
    
install_headers(my_headers,subdir: 'mincostflow')

pkg.generate(libmincostflow,description: 'min-cost-flow library')

libmincostflo_dep=declare_dependency(include_directories: my_includedir,
    link_with: libmincostflow)

subdir('test')
subdir('examples')

subdir('pymcf')
