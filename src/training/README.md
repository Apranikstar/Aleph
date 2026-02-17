```bash
 weaver --data-train '/eos/user/h/hfatehi/aleph/trainingFiles/aleph-temp2/Z*.root' \
 --data-test '/eos/user/h/hfatehi/aleph/trainingFiles/aleph-temp2/test/Z*.root' \
 --data-config /home/hfatehi/aleph/features.yaml \
 --network-config /home/hfatehi/aleph/model.py \
 --model-prefix /home/hfatehi/aleph/models/prefix \
 --gpus 0,1,2,3 \
 --batch-size 2048 \
 --start-lr 5e-3 \
 --num-epochs 20 \
 --optimizer ranger \
 --log /home/hfatehi/aleph/logs/train.log
 ```
