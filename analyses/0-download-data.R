# STEP 0 - Download TF and Cof data
# -----------------------------------------------------------------------------
#
# This script will download :
#   - Transcription factor and Co-factor data from the AnimalTFDB4 database
#
# The file will be stored in `data/derived-data/`.
#
# Ellis Franklin <ellis.franklin@univ-lorraine.fr>
# 2024/02/29


# Download the TF and Cof list from the AnimalTFDB4 database

dl_regulation_data(
  url_tf = "https://guolab.wchscu.cn/AnimalTFDB4_static/download/TF_list_final/Danio_rerio_TF",
  url_cof = "https://guolab.wchscu.cn/AnimalTFDB4_static/download/Cof_list_final/Danio_rerio_Cof",
  path = "data/derived-data/",
  filename_tf = "Danio_rerio_TF.txt",
  filename_cof = "Danio_rerio_Cof.txt",
  overwrite = FALSE
)

