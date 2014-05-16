FILE(REMOVE_RECURSE
  "CMakeFiles/clean-all"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/clean-all.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
