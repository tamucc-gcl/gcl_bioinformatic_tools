
install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
cmdstanr::install_cmdstan(cores = 20)
install.packages(c('shiny', 'shinyFiles', 'readr', 'readxl',
                   'writexl', 'dplyr', 'stringr', 'tidyr', 'DT',
                   'janitor', 'ggplot2', 'patchwork', 'jsonlite',
                   'yardstick', 'outliers', 'purrr', 'waiter',
                   'shinyBS', 'brms', 'callr', 'forcats', 'tidybayes',
                   'scales'))