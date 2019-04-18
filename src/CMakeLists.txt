target_sources(tersect
PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/tersect.c"
    "${CMAKE_CURRENT_LIST_DIR}/tersect_db.c"

    "${CMAKE_CURRENT_LIST_DIR}/build.c"
    "${CMAKE_CURRENT_LIST_DIR}/chroms.c"
    "${CMAKE_CURRENT_LIST_DIR}/distance.c"
    "${CMAKE_CURRENT_LIST_DIR}/rename.c"
    "${CMAKE_CURRENT_LIST_DIR}/samples.c"
    "${CMAKE_CURRENT_LIST_DIR}/view.c"

    "${CMAKE_CURRENT_LIST_DIR}/alleles.c"
    "${CMAKE_CURRENT_LIST_DIR}/ast.c"
    "${CMAKE_CURRENT_LIST_DIR}/bitarray.c"
    "${CMAKE_CURRENT_LIST_DIR}/errorc.c"
    "${CMAKE_CURRENT_LIST_DIR}/hashmap.c"
    "${CMAKE_CURRENT_LIST_DIR}/heap.c"
    "${CMAKE_CURRENT_LIST_DIR}/snv.c"
    "${CMAKE_CURRENT_LIST_DIR}/vcf_parser.c"
    "${CMAKE_CURRENT_LIST_DIR}/vcf_writer.c"

    "${FLEX_QueryScanner_OUTPUTS}"
    "${BISON_QueryParser_OUTPUTS}"

    "${VERSION_FILE}"
)
