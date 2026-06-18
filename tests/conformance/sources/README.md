# Maxima upstream (optional local clone)
#
#   git clone https://git.code.sourceforge.net/p/maxima/code maxima
#   cd maxima && git checkout master
#
# Extract:
#
#   python3 ../giac-rs/tests/conformance/scripts/extract_maxima_rtest.py \
#     -i maxima/tests/rtest_limit.mac \
#     -d limit \
#     -o ../giac-rs/tests/conformance/fixtures/maxima/limit.draft.json \
#     -p LIM-M \
#     --source-url 'https://sourceforge.net/p/maxima/code/ci/master/tree/tests/rtest_limit.mac'
#
# See scripts/manifests/maxima.yaml for batch jobs.
