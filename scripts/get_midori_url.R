## Get URL of most recent Midori Raw Database with all unique sequences for each species
## Outputs url to terminal
args <- commandArgs(trailingOnly = TRUE)
locus <- args[1]

Sys.setenv(TZ = "UTC")  # or "America/Chicago"

response <- httr::GET('https://www.reference-midori.info/download.php',
                      config = httr::config(ssl_verifypeer = FALSE))
httr::content(response, as = "text", encoding = "UTF-8") |>
  rvest::read_html() |>
  rvest::html_nodes('a') |>
  rvest::html_attr('href') |>
  stringr::str_subset('^/download/Databases') |>
  tibble::tibble(dl_url = _) |>
  dplyr::mutate(date = stringr::str_extract(dl_url, '20[0-9]{2}-[0-9]{2}-[0-9]{2}') |> lubridate::ymd()) |>
  dplyr::filter(date == max(date),
                stringr::str_detect(dl_url, locus),
                stringr::str_detect(dl_url, '/RAW/'),
                stringr::str_detect(dl_url, 'AA', negate = TRUE),
                stringr::str_detect(dl_url, 'uniq')) |>
  dplyr::mutate(dl_url = stringr::str_c('https://www.reference-midori.info', dl_url)) |>
  dplyr::pull(dl_url) |>
  cat()
