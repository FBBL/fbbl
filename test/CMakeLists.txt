cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

# custom function for testing
function(my_add_test name source dependency1 dependency2 expected_response)

  add_executable(${name} 
  	"${source}")

  target_link_libraries(${name} 
  	PRIVATE ${dependency1} ${dependency2}
  	)

  add_test(test_${name} ${name} 
  	WORKING_DIRECTORY ${OUTPUT_TEST_DIR}
  	)

  set_tests_properties(test_${name} PROPERTIES
    PASS_REGULAR_EXPRESSION "${expected_response}"
  	)

endfunction()

if(BUILD_FFT STREQUAL "ON")
  link_directories(${BINARY_DIR} ${FFTW3_BINARY_DIR})
else()
  link_directories(${BINARY_DIR})
endif()

if(BUILD_FFT STREQUAL "ON")
  # test_plain2_LF1_fft2_10_101_01
  my_add_test(plain2_LF1_fft2_10_101_01 "${TEST_DIR}/test_plain2_LF1_fft2_10_101_01.c" m fbbl "Test passed")
  # test_plain2_LF2_fft2_single_10_101_01
  my_add_test(plain2_LF2_fft2_single_10_101_01 "${TEST_DIR}/test_plain2_LF2_fft2_single_10_101_01.c" m fbbl "Test passed")
  # test_plain2_LF2_fft2_double_10_101_01
  my_add_test(plain2_LF2_fft2_double_10_101_01 "${TEST_DIR}/test_plain2_LF2_fft2_double_10_101_01.c" m fbbl "Test passed")
  # test_plain3_LF1_fft3_10_101_01
  my_add_test(plain3_LF1_fft3_10_101_01 "${TEST_DIR}/test_plain3_LF1_fft3_10_101_01.c" m fbbl "Test passed")
  # test_plain3_LF2_fft3_single_10_101_01
  my_add_test(plain3_LF2_fft3_single_10_101_01 "${TEST_DIR}/test_plain3_LF2_fft3_single_10_101_01.c" m fbbl "Test passed")
  # test_plain3_LF2_fft3_double_10_101_01
  my_add_test(plain3_LF2_fft3_double_10_101_01 "${TEST_DIR}/test_plain3_LF2_fft3_double_10_101_01.c" m fbbl "Test passed")
  # test_plain2_fft2_bruteforce_10_101_01
  my_add_test(plain2_fft2_bruteforce_10_101_01 "${TEST_DIR}/test_plain2_fft2_bruteforce_10_101_01.c" m fbbl "Test passed")
endif()

# test_smooth_lms_LF1_10_101_005
my_add_test(smooth_lms_LF1_10_101_005 "${TEST_DIR}/test_smooth_lms_LF1_fwht_10_101_005.c" m fbbl "Test passed")
# test_smooth_lms_unse_fwht_10_101_005
my_add_test(smooth_lms_unse_fwht_10_101_005 "${TEST_DIR}/test_smooth_lms_unse_fwht_10_101_005.c" m fbbl "Test passed")
# test_smooth_lms_fwht_bruteforce_10_101_01
my_add_test(smooth_lms_fwht_bruteforce_10_101_01 "${TEST_DIR}/test_smooth_lms_fwht_bruteforce_10_101_01.c" m fbbl "Test passed")
# test_smooth_lms_full_fwht_10_101_01
my_add_test(smooth_lms_full_fwht_10_101_01 "${TEST_DIR}/test_smooth_lms_full_fwht_10_101_01.c" m fbbl "Test passed")
# test_smooth_lms3_LF2_fwht_bruteforce_10_101_01
my_add_test(smooth_lms_LF2_10_101_005 "${TEST_DIR}/test_smooth_lms_LF2_fwht_10_101_005.c" m fbbl "Test passed")
# test_smooth_lms3_fwht_bruteforce_10_101_01
my_add_test(smooth_lms3_fwht_bruteforce_10_101_01 "${TEST_DIR}/test_smooth_lms3_fwht_bruteforce_10_101_01.c" m fbbl "Test passed")
# test_smooth_lms3_meta1_LF1_fwht_bruteforce_10_101_005
my_add_test(smooth_lms3_meta1_LF1_fwht_bruteforce_10_101_005 "${TEST_DIR}/test_smooth_lms3_meta1_LF1_fwht_bruteforce_10_101_005.c" m fbbl "Test passed")
# test_smooth_lms3_meta1_LF2_fwht_bruteforce_10_101_005
my_add_test(smooth_lms3_meta1_LF2_fwht_bruteforce_10_101_005 "${TEST_DIR}/test_smooth_lms3_meta1_LF2_fwht_bruteforce_10_101_005.c" m fbbl "Test passed")
# test_smooth_lms3_meta2_fwht_bruteforce_10_101_005
my_add_test(smooth_lms3_meta2_fwht_bruteforce_10_101_005 "${TEST_DIR}/test_smooth_lms3_meta2_fwht_bruteforce_10_101_005.c" m fbbl "Test passed")
# test_smooth_lms3_unse_meta1_fwht_bruteforce_10_101_005
my_add_test(smooth_lms3_unse_meta1_fwht_bruteforce_10_101_005 "${TEST_DIR}/test_smooth_lms3_unse_meta1_fwht_bruteforce_10_101_005.c" m fbbl "Test passed")
# test_smooth_lms_sample_ampli_fwht_10_101_01
my_add_test(smooth_lms_sample_ampli_fwht_10_101_01 "${TEST_DIR}/test_smooth_lms_sample_ampli_fwht_10_101_01.c" m fbbl "Test passed")
# test_smooth_lms_LF2_fwht_hybrid_10_101_005
my_add_test(smooth_lms_LF2_fwht_hybrid_10_101_005 "${TEST_DIR}/test_smooth_lms_LF2_fwht_hybrid_10_101_005.c" m fbbl "Test passed")
# test_lms_LF1_fwht_10_101_005
my_add_test(lms_LF1_fwht_10_101_005 "${TEST_DIR}/test_lms_LF1_fwht_10_101_005.c" m fbbl "Test passed")
# test_lms_LF2_fwht_10_101_005
my_add_test(lms_LF2_fwht_10_101_005 "${TEST_DIR}/test_lms_LF2_fwht_10_101_005.c" m fbbl "Test passed")
# test_lms3_LF2_fwht_bruteforce1_10_101_005
my_add_test(lms3_LF2_fwht_bruteforce1_10_101_005 "${TEST_DIR}/test_lms3_LF2_fwht_bruteforce1_10_101_005.c" m fbbl "Test passed")
# test_lms4_LF2_fwht_bruteforce2_10_101_005
my_add_test(lms4_LF2_fwht_bruteforce2_10_101_005 "${TEST_DIR}/test_lms4_LF2_fwht_bruteforce2_10_101_005.c" m fbbl "Test passed")
# test_lms5_LF2_fwht_bruteforce_10_101_005
my_add_test(lms5_LF2_fwht_bruteforce_10_101_005 "${TEST_DIR}/test_lms5_LF2_fwht_bruteforce_10_101_005.c" m fbbl "Test passed")
# test_lms6_LF2_fwht_bruteforce_10_101_005
my_add_test(lms6_LF2_fwht_bruteforce_10_101_005 "${TEST_DIR}/test_lms6_LF2_fwht_bruteforce_10_101_005.c" m fbbl "Test passed")
# test_coded2_LF1_fwht_10_101_01
#my_add_test(coded2_LF1_fwht_10_101_01 "${TEST_DIR}/test_coded2_LF1_fwht_10_101_01.c" m fbbl "Test passed")
# test_coded2_LF2_fwht_10_101_01
my_add_test(coded2_LF2_fwht_10_101_01 "${TEST_DIR}/test_coded2_LF2_fwht_10_101_01.c" m fbbl "Test passed")
# test_coded3_fwht_bruteforce_10_101_01
my_add_test(coded3_fwht_bruteforce_10_101_01 "${TEST_DIR}/test_coded3_fwht_bruteforce_10_101_01.c" m fbbl "Test passed")
# test_coded4_fwht_bruteforce_10_101_01
my_add_test(coded4_fwht_bruteforce_10_101_01 "${TEST_DIR}/test_coded4_fwht_bruteforce_10_101_01.c" m fbbl "Test passed")
# test_coded42concat_fwht_bruteforce_10_101_01
my_add_test(coded42concat_fwht_bruteforce_10_101_01 "${TEST_DIR}/test_coded42concat_fwht_bruteforce_10_101_01.c" m fbbl "Test passed")
# test_utils
my_add_test(utils "${TEST_DIR}/test_utils.c" m fbbl "Test passed")

