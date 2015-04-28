# message(STATUS "Making test data")
# add_custom_command(TARGET test.vcf.gz
#                    POST_BUILD
#                    COMMAND "cat test.vcf | bgzip > test.vcf.gz")

# add_custom_command(TARGET test.vcf.gz.tbi
#                    POST_BUILD
#                    COMMAND "tabix -p vcf -f test.vcf.gz")