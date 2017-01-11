import sys
import os
import logging
logging.getLogger().setLevel(logging.INFO)

scriptDir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(os.path.abspath(os.path.join(scriptDir, '..', 'python', 'Tools')))

from fastasize import fastaContigLengths, calculateLength

def main():
    locations = "chrMT chrY:1-10"
    fastacontiglengths = "{'chrY': 59373566, 'chrM': 16571}"
    total_length = calculateLength(fastacontiglengths, locations)
    if total_length == 10:
        logging.info("fastasize test SUCCEEDED!")
    else:
        logging.error("fastasize test FAILED!")


if __name__ == "__main__":
    main()
