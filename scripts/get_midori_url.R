## Get URL of most recent Midori Raw Database with all unique sequences for each species
## Outputs url to terminal
rvest::read_html('https://www.reference-midori.info/download.php') |>
  rvest::html_nodes('a') |>
  rvest::html_attr('href') |>
  stringr::str_subset('^/download/Databases') |>
  tibble::tibble(dl_url = _) |>
  dplyr::mutate(date = stringr::str_extract(dl_url, '20[0-9]{2}-[0-9]{2}-[0-9]{2}') |> lubridate::ymd()) |>
  dplyr::filter(date == max(date),
                stringr::str_detect(dl_url, 'CO1'),
                stringr::str_detect(dl_url, '/RAW/'),
                stringr::str_detect(dl_url, 'AA', negate = TRUE),
                stringr::str_detect(dl_url, 'uniq')) |>
  dplyr::mutate(dl_url = stringr::str_c('https://www.reference-midori.info', dl_url)) |>
  dplyr::pull(dl_url) |>
  cat()
