	B`??",@B`??",@!B`??",@	c??%?q??c??%?q??!c??%?q??"{
=type.googleapis.com/tensorflow.profiler.PerGenericStepDetails:B`??",@bX9?ȶ?Aj?t??+@Y#??~j???rEagerKernelExecute 0*	      [@2l
5Iterator::Model::ParallelMapV2::Zip[1]::ForeverRepeat???S㥫?!      I@)?A`??"??1??^B{?H@:Preprocessing2F
Iterator::ModelˡE?????!???^B?B@)???S㥛?1      9@:Preprocessing2U
Iterator::Model::ParallelMapV2y?&1???!?^B{	?)@)y?&1???1?^B{	?)@:Preprocessing2v
?Iterator::Model::ParallelMapV2::Zip[0]::FlatMap[0]::Concatenate?I+???!?%???^$@)?~j?t?x?1??8??8@:Preprocessing2?
OIterator::Model::ParallelMapV2::Zip[0]::FlatMap[0]::Concatenate[0]::TensorSlice{?G?zt?!Lh/???@){?G?zt?1Lh/???@:Preprocessing2f
/Iterator::Model::ParallelMapV2::Zip[0]::FlatMap9??v????!/????(@)????Mb`?1?Kh/???:Preprocessing2x
AIterator::Model::ParallelMapV2::Zip[1]::ForeverRepeat::FromTensor????MbP?!?Kh/???)????MbP?1?Kh/???:Preprocessing:?
]Enqueuing data: you may want to combine small input data chunks into fewer but larger chunks.
?Data preprocessing: you may increase num_parallel_calls in <a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#map" target="_blank">Dataset map()</a> or preprocess the data OFFLINE.
?Reading data from files in advance: you may tune parameters in the following tf.data API (<a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#prefetch" target="_blank">prefetch size</a>, <a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#interleave" target="_blank">interleave cycle_length</a>, <a href="https://www.tensorflow.org/api_docs/python/tf/data/TFRecordDataset#class_tfrecorddataset" target="_blank">reader buffer_size</a>)
?Reading data from files on demand: you should read data IN ADVANCE using the following tf.data API (<a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#prefetch" target="_blank">prefetch</a>, <a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#interleave" target="_blank">interleave</a>, <a href="https://www.tensorflow.org/api_docs/python/tf/data/TFRecordDataset#class_tfrecorddataset" target="_blank">reader buffer</a>)
?Other data reading or processing: you may consider using the <a href="https://www.tensorflow.org/programmers_guide/datasets" target="_blank">tf.data API</a> (if you are not using it now)?
:type.googleapis.com/tensorflow.profiler.BottleneckAnalysis?
device?Your program is NOT input-bound because only 0.6% of the total step time sampled is waiting for input. Therefore, you should focus on reducing other time.no*no9b??%?q??I???o?X@Zno#You may skip the rest of this page.B?
@type.googleapis.com/tensorflow.profiler.GenericStepTimeBreakdown?
	bX9?ȶ?bX9?ȶ?!bX9?ȶ?      ??!       "      ??!       *      ??!       2	j?t??+@j?t??+@!j?t??+@:      ??!       B      ??!       J	#??~j???#??~j???!#??~j???R      ??!       Z	#??~j???#??~j???!#??~j???b      ??!       JCPU_ONLYYb??%?q??b q???o?X@