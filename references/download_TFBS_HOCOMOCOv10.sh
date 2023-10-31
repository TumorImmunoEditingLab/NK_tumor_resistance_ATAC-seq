#!/bin/bash

filename="TFBS_mm10_PWMScan_HOCOMOCOv10.tar.gz"
curl -0 https://www.embl.de/download/zaugg/diffTF/TFBS/$filename -o $filename  && tar xvzf $filename  && rm $filename