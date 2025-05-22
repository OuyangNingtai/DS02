<a name="br1"></a> 

Fast Approximations and Coresets for

(k, ℓ)-Median under Dynamic Time Warping

Jacobus Conradi [#](mailto:jacobus.conradi@gmx.de)

University of Bonn, Bonn, Germany

Benedikt Kolbe [#](mailto:benedikt.kolbe@physik.hu-berlin.de)

Hausdorﬀ Center for Mathematics, University of Bonn, Bonn, Germany

Ioannis Psarros [#](mailto:ipsarros@athenarc.gr)

Athena Research Center, Marousi, Greece

Dennis Rohde

University of Bonn, Bonn, Germany

Abstract

We present algorithms for the computation of ε-coresets for k-median clustering of point sequences

in R<sub>d</sub> under the p-dynamic time warping (DTW) distance. Coresets under DTW have not been

investigated before, and the analysis is not directly accessible to existing methods as DTW is not a

metric. The three main ingredients that allow our construction of coresets are the adaptation of

the ε-coreset framework of sensitivity sampling, bounds on the VC dimension of approximations

to the range spaces of balls under DTW, and new approximation algorithms for the k-median

problem under DTW. We achieve our results by investigating approximations of DTW that provide

a trade-oﬀ between the provided accuracy and amenability to known techniques. In particular, we

observe that given n curves under DTW, one can directly construct a metric that approximates

DTW on this set, permitting the use of the wealth of results on metric spaces for clustering purposes.

The resulting approximations are the ﬁrst with polynomial running time and achieve a very similar

approximation factor as state-of-the-art techniques. We apply our results to produce a practical

algorithm approximating (k, ℓ)-median clustering under DTW.

2012 ACM Subject Classiﬁcation Theory of computation → Design and analysis of algorithms

Keywords and phrases Dynamic time warping, coreset, median clustering, approximation algorithm

Funding Jacobus Conradi: Partially funded by the Deutsche Forschungsgemeinschaft (DFG, German

Research Foundation) - 313421352 (FOR 2535 Anticipating Human Behavior) and the iBehave

Network: Sponsored by the Ministry of Culture and Science of the State of North Rhine-Westphalia.

Aﬃliated with Lamarr Institute for Machine Learning and Artiﬁcial Intelligence.

Benedikt Kolbe: Partially funded by the Deutsche Forschungsgemeinschaft (DFG, German Research

Foundation) – 459420781.

Ioannis Psarros: Supported by project MIS 5154714 of the National Recovery and Resilience Plan

Greece 2.0 funded by the European Union under the NextGenerationEU Program.

Dennis Rohde: Partially supported by the Hausdorﬀ Center for Mathematics. Partially funded

by the iBehave Network: Sponsored by the Ministry of Culture and Science of the State of North

Rhine-Westphalia.

Acknowledgements We thank Anne Driemel of the University of Bonn for detailed discussion and

guidance.

1

Introduction

One of the core challenges of contemporary data analysis is the handling of massive data

sets. A powerful approach to clustering problems involving such sets is data reduction, and

ε-coresets oﬀer a popular approach that has received substantial attention. An ε-coreset is a

problem-speciﬁc condensate of the given input set of reduced size which captures its core

![ref1]![ref1]![ref1]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAYAhQDASIAAhEBAxEB/8QAHQABAQABBQEBAAAAAAAAAAAAAAgBAwQFBgcCCv/EADIQAAECAgoCAgECBgMAAAAAAAEAAgUGAwgREhhSWJKX1gQhBzFBFCITQlFhcYElKDf/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A/fwimuaaxHmy5MPny/4/wNWDmWj8LyDRCZJUkqW4nLfkgOIs8bzfOnWF+S8my03vCZYCD7J9cfiZiJt/61VoRYbP/P5V92fkWT+bQfwgqRFLeJmI3mA1bqzbLzg22lkGV2A22+mls+Ptef5WkAGw2uC9F+PPljyp+/XjyPi75U+PneFaWP8AkiX4RAqLyx+P0phkwx11J/e8xn2P6+g9fRYHsA/1AWUBERAREQEXxSEhpIJBHv0A4n+wBIBP+11GZ5pppZhMRjLIFH5m/QUTXsgMqeF4sRmHzXG8bvi+H5vmwzx3v/bZdd5jWlzgLwQdxRS3iYiOmutBx/KvfkxMRHTXWg4/lXvyCpEUt4mIjprrQcfyr35MTER011oOP5V78gqRFLeJiI6a60HH8q9+TExEdNdaDj+Ve/IKkRS3iYiOmutBx/KvfkxMRHTXWg4/lXvyCpEUt4mIjprrQcfyr35MTER011oOP5V78gqRFLeJiI6a60HH8q9+TExEdNdaDj+Ve/IKkRS3iYiOmutBx/KvfkxMRHTXWg4/lXvyCpEUt4mIjprrQcfyr35MTER011oOP5V78gqRFLeJiI6a60HH8q9+TExEdNdaDj+Ve/IKkRS3iYiOmutBx/KvfkxMRHTXWg4/lXvyCpEUt4mIjprrQcfyr35MTER011oOP5V78gqRFLeJiI6a60HH8q9+TExEdNdaDj+Ve/IKkRS3iYiOmutBx/KvfkxMRHTXWg4/lXvyCpEUt4mIjprrQcfyr35MTER011oOP5V78gqRFLeJiI6a60HH8q9+TExEdNdaDj+Ve/IKkRS3iYiOmutBx/KvfkxMRHTXWg4/lXvyCpEUt4mIjprrQcfyr35MTER011oOP5V78gqRFLeJiI6a60HH8q9+TExEdNdaDj+Ve/IKkRS3iYiOmutBx/KvfkxMRHTXWg4/lXvyCpEUt4mIjprrQcfyr35MTER011oOP5V78gqRFLeJiI6a60HH8q9+TExEdNdaDj+Ve/IKkRS3iYiOmutBx/KvfkxMRHTXWg4/lXvyCpEUt4mIjprrQcfyr35MTER011oOP5V78gqRFLeJiI6a60HH8q9+TExEdNdaDj+Ve/IKkRS3iYiOmutBx/KvfkxMRHTXWg4/lXvyCpEUt4mIjprrQcfyr35MTER011oOP5V78gqRFLeJiI6a60HH8q9+WlS1nImxrntq1VoKUsu20NH8fSq6mdeBs9O+QGNaPX3ePr8IKpReVyr8lRCZYN48X8n44n+V6XyHP/4eZoVD/CjHjtbZdPleN4sW8+goy8H0GeVSewffr2QegU/hP8tjqLyRRUrBSNez+IwUrHXbbHGheQ1jhb6IJstNi3tFQNowQ11JY5xfYXkhpdZ+1g+msFn7Wj0ERBl9EXNIbSPYT/NaXev6WEhbUeITS/xSyhD6M/se1jS5/wBe3EgFh/wXfn/ZEHIC2wW/dnv/ACiIgIiICIiDR8ijNLRlg9EkewS0t+/bXD2HD8EfS2NJDmucykLKPyHUNFco/wBS0UlIX+vbvIIdSD6ykevyiIOQDPQtc8GwWgPNgWbgzP3lEQLgzP3lLgzP3lEQLgzP3lLgzP3lEQLgzP3lLgzP3lEQLgzP3lLgzP3lEQLgzP3lLgzP3lEQLgzP3lLgzP3lEQLgzP3lLgzP3lEQLgzP3lLgzP3lEQLgzP3lLgzP3lEQLgzP3lLgzP3lEQLgzP3lLgzP3lEQLgzP3lLgzP3lEQLgzP3lLgzP3lEQLgzP3lLgzP3lEQLgzP3lLgzP3lEQLgzP3lLgzP3lEQLgzP3lLgzP3lEQLgzP3lLgzP3lEQLgzP3lLgzP3lEQLgzP3lLgzP3lEQLgzP3lLgzP3lEQLgzP3lLgzP3lEQLgzP3lLgzP3lEQLgzP3lLgzP3lEQLgzP3lLgzP3lEQLgzP3laLqKkDnOY5zw4NBbSUjrjQLbS1vv2fs+hb+SiINajaWNuucXm0+yLPv8WBERB//9k=)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAWABkDASIAAhEBAxEB/8QAGwAAAQQDAAAAAAAAAAAAAAAAAAYHCAkCBQr/xAAzEAABAgUDAQYBDQAAAAAAAAABAhEDBAUGIQAHEjEIExZBUWHRFCMyM0JSU3GRoaKx4f/EABkBAAIDAQAAAAAAAAAAAAAAAAYIBQcJCv/EACURAAEEAQMEAgMAAAAAAAAAAAIBAwQFBgcRIQASMUEIIhMUUf/aAAwDAQACEQMRAD8A6itzdzr/AKbuFdUlJ3fXqbTKdWIiZWTQsKSpKjFBQg96kBA4hsDqceWkTC3Z3JCAvxvcEMxXi8QtvrC+fnnyXf0/LOsd2wIe5l5fa7urz3ErdbKhxEBH0n6c1OSQ+D5Bm8KeKlhwwJ92z0f0YdPIAa5YtT9TtQo2oOXomYZWIjlk5BFLEdkROOPttwi+NuOeOUTp/MVxXGX8Zr3nq+M465Gjmbhx2iNTJoCLciBfaqvKr59+3coW5+4M7cFvy0xe1XjS0eu0eFMwZlPeJjQIlRlkqhJaOSgqcAqbCeQbJ1bPzh/iK/l8NUrWsCbotsl3FwUUukkAkVKVGfVurH39H1dE6/f9P81oL8B87vbXEM7k2k+daSUvqgP2bF788hQSsJRDvRV+o+BT1x/V6pXXKrra2wx1mtjtRmSr5ZkLLbbXeavtIpH2CiEW2ybrvsnG/HUML47K9fua7a9cEvddHloFXnJqahy8eSnVxoKZhaFJStaCEKUnhnjglmZtJY9ju5ipR8ZUPJJzIT5OfXIGQA+G/vRo0X5l8cNFLLIrqbNwWG/Jk2smS+6trkAK4+a/ZxRbtgBFXfwIiP8AB6jKvUHL4dZGixrgm47bLQg3+lWmgiIAgohHDI+ERPJKq7c9bKidkW46dWaTUIl30SKiQqchOqQmRnwuIiUmoMwpCSVcQpYhlIJDAlyWfU+vkJ++P3+GjRphPjvpRp9hNJfw8WxxipjTLGG/Jabm2khHXW434wNSmTpJiogqpsBCK+VRV56B82yS7vpUB22nFLOPGJplVYjM9gE42RCgx2WRXdfZIq+kXbjr/9k=)![](data:image/jpeg;base64,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)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAbABUDASIAAhEBAxEB/8QAGgAAAgIDAAAAAAAAAAAAAAAAAAgGBwUJCv/EACkQAAEFAAEDAwIHAAAAAAAAAAMBAgQFBgcREhMACBQhMhYXIiMzQVL/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A6+eVOXPczM5+/LjgXH4O0z3HNFk9VyLH2WkNUzNLA17rocOmoTDqZiZ6XXLRS3SLZ47Ty/KC34Q/H1fRtBuPe7z3zbrMLFJjfbtjuJHUM+Wh4cnk631lwcE8dMYUtS4pmhw9zFZZzJP7cRQSI1X/ACoTq3YfGweSqd3p+SCk8V1vqnL0Vy+TZiWvfFxrrh9M6LFIISoqrf2CzCoZUd2gTtTp9ZrGVxCFM00QwzPaBjxKN4pLY6Pa1rjjVehS9UexF7vjIxzE7+/q0Kv4c5B0F7P2nHWyYEvIXGBM9H1FnDF4Ki/haaHNm5y/gCc8jq49vCr5EywzyklpnzEZASzs+nyVPUb4qfh7HlnnvT5DkCu11rPscRm9lm6NjTycLe4utuqktfcy2yFU86xIWSVWujRlC6GVvQncrmnoF20eE4x5X96GvwnNXh1cHP8AEmKl8XYPRSSuyk0lzLvG8g2UKhKN0KROivr8mx7myWnGw7URr0e5W3PM4b1XEhx2Xt8s6ikzcYL7fScQy6ZEzd9MA5iKHJ3A5jncdTbUZZMi1kRKe+W6lxoKFZH8CPVidDk81optGS8o620JUaSJfVhJsYZnwbmtGcdfZRnuTuHLiMlyWgK1eo0ORE+71nJLUYEhWJ2k+aYXen3eN6uVzOv+XK1vVP76J6BHPa3b8e23N/uu0OOpPw/XX1lxLYXdNJqhVVvXbguf0xNpE0UMRjMHpQ3BXjvEaY7EnNIjJB0/Wp6ayxz1Hm9Fb39DVQqi61QKpmjs4IGAmXTaEBYdMliZiI+SlZFlSARPIq+EZiNb9HL6PQf/2Q==)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAKABUDASIAAhEBAxEB/8QAGQAAAgMBAAAAAAAAAAAAAAAAAAYHCAkK/8QAKhAAAQMCBQMCBwAAAAAAAAAAAQIDBQQGAAcREiEIEzEVUQkUFkFScaH/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A6180s9us2BzYumyssemKBzGtyFZhpKguT66NsNuRMz6n2IovXCxFUC7kpfTwZFmCq5KhbQ633q1JUyFq6eoj4jS2ULV0M20gLaW4sqzutTuMPd5tKKcvIlFIJU2pwqKFqQko2AnUY0eryfl5hWp3JpqFSTrylW9Y1SfIOhI1Gh0JH3wwLSnvoRtTt/HQbfHt4/mArf0+X1n7fELcNXnlktG5TTdDLU9NDQtDdMbc7VfGuUaXna8yFA++0FKqCQWyvcCo6jjBifiTtHJ4cfA5PgPuAD9AcD2wYD//2Q==)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAOABUDASIAAhEBAxEB/8QAGQAAAgMBAAAAAAAAAAAAAAAAAAcBBQgK/8QAIRAAAgMBAQACAwEBAAAAAAAAAwQCBQYBBwgUABIVERP/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A66vbfb/k5Y+w3Xh3xu85zLNpl6DLbDReh7K7Y5mRp25bQMMh/OXqCzWsn4qFOG1+6aI4pk59En7d7F0ZH3K8/t1uP9h8+Y831loZlapskbSWj8usbTkZmVzFZ6AetomG9K1WDdsjVZ8ymspCsfFB9qYh9NFtgdfkPdHvWKHRVjtD6OLDYTW4q2SIL6qWfDpWEdFS3QING/oJysDr9oSqLV9oNuTDNoqSuXEy69dmKrW5ixymnVFcUdsqNewXLzg/+hQtrsLGCP8AUkAEWaAF1UnJSksyuCQ+d7Hkohn/AOPXoXqd/t/kNjPYjZX7fnHoVRXZCecUMqOWM0uaT1FBCxnMhetWS9XYpAbYjEUSHgSURx5L/OH5ZfHzzHU4Db+9P6LWB1PdTqcUSnckoQVitSZvA0+aqlLghO947ZhRr1htWXCknYlhN00QmPMUT8D/2Q==)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAGAAYDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAcEAABBQEBAQAAAAAAAAAAAAAEAQIDBQYABxH/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8At7xeC9DrN/6nob/1a40ue1V4EXlctMCEALjggK8evkBEmBghkexz4HIiNcjCm/Dz0ntZyipXOcH/2Q==)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAbABUDASIAAhEBAxEB/8QAGgABAQACAwAAAAAAAAAAAAAACAACBQYJCv/EACQQAAEFAQACAgIDAQAAAAAAAAQBAgMFBgcACBESExQVISIk/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/APZn1LQ+0l5v34jkORx+dxcK5axtOv646e/+g1o61bbh5zAMDrorqwqEFFknaRp6tqoXA1Hf7VWl/ofTvfnlHZ8niqgzmHeM7YZa2201bXZKfn+j0I2TNpK7W1oCsudRDVFpLpqsil+HmsMZEYyVwyxtWRJexNF3AAF+34f1zJ4jTmy1NbdVnUXwmc5cDChboJ6pbNsrKG7YssyFSBALPfMdB+++P+NGRxg9attmC+0aG+7p7MYve9/zOeu8lU5LHHlUPN8Tk0Pq10DqBbEepoTjrk0Wk/lyoVhk/wCAX9RTo3zOgDsXxW1oOkVUl1ST6ICUQmWruKcpI63R0duGqxmU+lASUhg9iBKjoHJDOQOr2zfgImYn28vM8ti8eFqNb0utk/PcdJr8jDcEw2DrWiJDyIdkLRkUccayCwQywW5TpyB1RT1WGWX+4m+XgDDouL4JofYcmT2V6GyU8SOsuec8uvNFp6/jw9Zn5C/xXulqrt4eA1mvtVPiZcVtxFaQhNDHQD7NLLReVbPn/oLqMda0V5VetIdVZirV2p9AzmGfu4x5GuIJYBc5/wDUtKyZVFjVktKVGSis+g6q1zk8UWhxOM3B88OyyOX1MYyxwDJf0FTauHjPZO0tIXmiTPifMkEXzKxzZWLG10b2OT58Jub9CfUbG63NaSg4vSw2oF1G8ZbG92N/XfZ0BKqpFHf6OzpDPhWp9UMrp0b/AGjURFX5Dc+odxy5tVqcPxOy6Kdx3GB4x2FJ1q6J1Qo16BbTFw5DSbZ0Wiv6qFABGyQOkmpqZFFhpfxjlSo68bFPT1NfXBV9fV1wNeCGMKDXhBDCgAjQsVkQwQUEUYwg8bERrIR4o4mta1rWoiIiXgf/2Q==)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAaABUDASIAAhEBAxEB/8QAGgAAAQUBAAAAAAAAAAAAAAAAAAQGBwgJCv/EACUQAAEFAQACAgMAAwEAAAAAAAQBAgMFBgcIEQATEhQVMUFRkf/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDr16n1byhse/z854Lm+aEZ3mufyGq3g2801hUk6wbYOu4xaSpICzVwmZlrf4RLpbhqWchH7cbUBZ9Xt0N5fT+cfde16vGWJ2P8b8Tyv+Ee/wDnV5XVrHY280J0dMUPYWUOAivcXaCtsibAYlQCgyoqpFGm/NzorUbvpHjrxvrrStvs6LLdP6+Fm86SPaGWBNcSJknWjqFx6RwEV+ehV17ZfrHXTq0aw9Sp98qiL9b26T5GcV5Dmq3a9F6JR1GW0csImcMEnI0a36ktV48tLXUEVrY28HpqJKfXiEgxOkhZ+yjp4UkBXx3f6G6tN1zfZJCTveWzZ6DQ3QcbYKjR1unEPOzV8LB7e+ssLICumKus6ilj54uVleLcXMafuKfGlwy65dvtL1Pq3NNxJpF15+ZoNLUy1tpTn5i1wYNnTIJY01uDX3wJJakkPc+1rx3lOHfJA6ZrZHNPgRf5DdQ55lr7aB6jxf2/XNDVZGlJoL6m41Z72u2BVk22Z/CS9qs5YQV8lJ9DHTzGlxxjNsmuFdA50jnZ38m407m6Ynyx8wsRf1vOclmVruY8xpBLLYZXiVXG8J2Zi0uSVtzqJHwgQERCWhRDgQ3RTRbNx1mTRywb3zqqW1e1FVGrGeit/wBKiOF9Iqf4VE9r/wCr/wB+Kzwgyx2BlCDEiTSsjmFIgimGlj/F6/XLBIx0UjPaIv4uarfae/XwK1ePt7zrdu2XYeZ4ewz1P05+bOn2FxWlZYnoslUCdFFdJmD4giw2DNNkRLh4MMWj/dQ+GctjftQ+WWrxhhhkGGHgHHge+KCCCKOKGGJi/iyOKKNrWRxtRERrGNa1qJ6RET4fA//Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAoABUDASIAAhEBAxEB/8QAGgAAAwADAQAAAAAAAAAAAAAAAAcIAgMJCv/EACUQAAEFAAEEAgMBAQAAAAAAAAUBAgMEBgcIERITABQVFiEXIv/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwD1yc488dSl3mE3wP028bgr2gzudzGwO8g6wvKuVqUispRjMxLQiGqsBEglCb031vO9aV5F+nJ/ezpyPOuhjNCMhzJx3c451Ji2+iJJjSDtLxlZLPajh2ai38lAL5astAy9ajAKDaleAbaX7syt7fNVzjjbZfn7S8yZ08KKZjeBMXjtbii49w2cILxbz80mxF6Vtq2ltleM5L7QCiIPySeDvylZK/aRrPlwnL+EextkdrMRrR/ujsUnx2R1yjZVklK3TuM8kRiKz21rDWeauRHeLeyooKHhTkLkwtzH1H4Dk92ZiE4EzibnHCha8kVn9O2FbT3KMRmysj0ukIK4qm2SZsUDfN0ipH/32afM+D+KtbxdyTzSd12nraQFsKvGY3EFbr1TRoHxI3S0JB+mR7pUtlR7StOKY218P5qR8tlaFHwSNx8CiLUdL03a0roGQQj5kuVnR+iP6ssb/W2WXyk8GxRMs+a+Dkc16/xqJ/eUmv536Z8mOxPHnSftz6bOAgefhchwKHbtMGzU221HtfyBiYSAhSUcyV5Ihw9DFPxa69Mk7liWN9S9Q/WHwR08mRmO5CI6OY5o6EkrqGNDu0hmtRII2Bbt932KKDUuor3VLXawk/1bDvWz0dnqyHd9FvRRicYfyNKBf9RQQoK7kxdfScg8lV2QTLDqCD2vHyGI4HXY2E7XjVRtgnVRO/t7IFI9N+25z0vHYYr1G4oJx9r7YUFZ+uOPSWpyVyxTkeamIAZBVRcnYr2frt/C/fMLAs0kS3n+jykPjE4t5NAb6sYrwCtFlS+empVy2R2oP9d0guK6yw8aRsCXWrvpoG4qticbKlh3uirSL2TxTufAUPNXKuywCaF4DgLW8pWHZ2abBXMnSYThXR02vS+D2L3U3OxleF849RZFqmVLo6+36lL6KOnjnp1xnOOAvGOpjqT4sL77d8sDmjkzeKzyl9rxSFqzMshMtRxE8rEiFXmWLLjqMK1f160IFV0YS/IJJUPh8DoTxZotjvTm50J7FksliZH56rhh2vEqF2jrdauTbrpjQ1ZJ/TQfbUWoORLMvugW0vjH4ojj4fD4H//Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAUABoDASIAAhEBAxEB/8QAGgAAAgIDAAAAAAAAAAAAAAAAAAgHCQIGCv/EACgQAAEFAAEEAQMFAQAAAAAAAAQBAgMFBgcACBITERQVIQkWIiMxQf/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDvMuzX1VfZWcCOIkrqo4weJy+MU7oBpJ3Dvn/KseR6WNjf4KkSo5VRfL8Ll2o7flnk7jWt5U5PKqBoORRoNTjcPUCpG/CYy3b9XTVV7dK9v7m0aDSitNL+kChFlgIZFISk/kkrYPfVPIjuQoBq0qITG76847txj4mqKURQMiac8SGPzYSPM01jVarHJMiIjkX4ROk55E/T/fstSbosD3I828OZew9D6vjzEXloLkqIfxc4gSlqZrIeCjr55EiePUAiCgVTImDAQQxNRiAz3DXIt5puQucePrmzqdDHxdo6Aer0tX4iusQtjVTaP7PY1MLiRwzsgrmZ55jLAiW2UeQ8katmkcHGxfSxcacY4TtQ4j0rs6HYW32oC62Wy0R73G63a3I8BVra3FzZGyynWJc0qluhaUTKgzJlhFRsKI3qdcdpY9lkcrrwYfQFq85R6QOCfybPCLe1gtoPFM1EXxljiKYyRP8Aj0VOgqp5a7J8hquUN3ro+Y+4XJyavV2OnsKDGcgVNFmx7izeimEAgJkiSI/b6o0VCDiXfDG/yVU+etAf2FZBZFHZz93VRf0yvWVnLo/tVzJ4PFV88s+NP9VF8Y0+UVU6OjoMZOwbJHxLGX3A91E8MkpohA83KdRMOWHM6UeYQmGbFSMlglhc6N7fhHK1VRHp1apjqKHJ5HK5arNsnVmazdHn65xZEcpTgaasFrhHEysghZIQo40azSMiiY+Tyc2NiKjUOjoP/9k=)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAXAA8DASIAAhEBAxEB/8QAGAAAAwEBAAAAAAAAAAAAAAAAAAYHCAr/xAAhEAACAgMBAAIDAQAAAAAAAAAEBQMGAQIHCAAUEhMVFv/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDsm7J1W3cZb9P6BbrdV9uWVLny+et0NUpCPvj21z7SwHOSMktl+wwqovZdosTx6F6uhj2E8hi/dZFCbEOP9s9Hr03BnPabBSLCw77s6+jWElf0q5NOWiU1/eVFl3O3NONZsGa9HEpaVndUAGlncTzaPj91kMBz9b6pz3t/pr/ONuE1q1r+Yph3tw61OGvz9m6wzxwVOkjF42/J9EpE3sZNlTt8foUlxpN9YsyzaZ1Xbog4bD615yf14aYPt59fvh3KGCuw2smnxVYAQ1OUpASTY0Vh2RhVSCXLmOALVbvKKeTobMdoLiYLZ1hL6Ep+lmtnAyOeu4z5MEs6FeV264bZ1BJjSdygbV6dDn+rYsE7SWIi0mtIt9wV387QPGC8T58575l9A9d7lUfS3pt1WKczpVcsKCj8roOBmYdc0ciTJSiS3JezqAyFyvIIbzfoOycKx3hEilgXfZCkPh8D/9k=)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAGAAYDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAcEAABBQEBAQAAAAAAAAAAAAAEAQIDBQYABxH/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8At7xeC9DrN/6nob/1a40ue1V4EXlctMCEALjggK8evkBEmBghkexz4HIiNcjCm/Dz0ntZyipXOcH/2Q==)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAmABUDASIAAhEBAxEB/8QAGQAAAwEBAQAAAAAAAAAAAAAAAAcIBgQK/8QAJxAAAQQCAQUAAQUBAAAAAAAAAwECBAUGBxEACBITFBUJFiEiI0H/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A9h6793IXuRtNCRNFQ52L0Y66/udqEz0gK6uxm7LKHWvsK9cUIw17YpXy3x6b8kNsr5Sp9gvWiu4N9dyXcJq7OZGN677R883LjkWkrZiZhjdyWuhTZsxpVnVo4LKCzT21yjE1jvrd7GncvgPx4cwKDbWLVe/tlaMyCcUFulTj2yKiXkd82Ya9hZUW3ZcUlHBNAgigUeGLVQVjoljNe3847/ASJyShZthU1cENpJsI4YSDc4lsSYEUR4Vb7FOh3PUZTuYNXgYxy+bGlVF/joJ17b99Xu5YmXwdgaxyTS+cYXLqA2+H5Qx8uIGBeBsD0xqixcCvWY5wa+T+ST4xIAqxkRX+fKHW81DuXHtz3edysMjR7nX+PpjETHtkQJKSaLM7WZGtX5NWVqLGG8UnEJUWHDs0cY6OLZB8fWjVRx0CG7nML7V9/Zbi+jtvZVVxNnJIh3WK1kGyNEzWsfaewUaPWSBwiIyqvSR0+qO84GSVrgcqqD5ScIP6YeUT63DsF2D3U7BzzS+K2dZJTVsipfW0kSDUCkx69lNL/cs9lTNhhlkFDmjgHVGEKnqai8dXphCY4u5dp2Elauv2fLqcUg3WP1GTyrcpMGqH3y4LdGpZVbWChGsnTcgSUyOaYxvzjR514by63B9YP8U4jxGEHHitM+JBY1XMR0aUETCtVQKjUC9Gv58iLwz/AKCl0fr/AFrqvX8TXWnwGTFsOtbWgltjSVNYrkdeUYbsuQzVAB0++MdGOsJThIpCorvJ/P8AB1lsAyO/n9wPcTjqXtjY1GNV2nlra2XOeyrpz21RlxbFK1onyk+uaSIB9u54YpHmBGV7SLwrDoMf3Kdv8nJ7OBvfTt3F173F4RTym47l0kZ20mVUY1C8mF7DFAYaXa44YgxOHzHsPiRZPqgm+l/EihvP1N5i/Ofb3bmKGL5xzzRcWv2S3RpPmP3wB/tto4swTUe5iIdGefivuTx56OjoL20Poip7fcbngj3c7L84zSYO72NsO7D55Bm+QBQvFnaGJJkHd6nTZfzDLKkKFkgjWv8A7LydHR0H/9k=)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAcABoDASIAAhEBAxEB/8QAGAAAAwEBAAAAAAAAAAAAAAAABQcICQr/xAAnEAABBQABBAEEAwEAAAAAAAADAQIEBQYHAAgSEwkRISJBFBUxFv/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDqGPwf3P8AOG/5o02Y7w+TOJqCi5k1uUosVBqDz6oFFUrEdWmqV/vq1HDOkgqM8Q+LkGn0cv36Lg7Ke8EgmEf8h3KI3ETzVi5k8Z7fL7+JRf8AXF8Sp/j/AM1+qp+uq11dlzrb8k3WNxNVn+P8LUQc/pj8tXkZNRIureWWe20zNNhnJVBe1BBA+bo36Fha57orGVctJTnBQ/OFn8hVhy9qqztzBxNU8X09bROgWXIUQnuuLmbDKa9WFMCGQQqVkoXpPHeFiDcVjWkcn36BP7btI7v8jk9Hq0+QPla0fmqO3ux1QaSQN9matrpMqPB8k0xV+kg4hjVEG/6+X+L1opxJBurnirjK40FvNkX1tx9jLK7kS2qyUe3nZytlWRpLHOVzZBZhTPM1yq5pHORVVU6VXatu+etPG5AwfclR5QG+wU2hETRYk5jZbWVWrqX3cUsKNIr4axC1kYgK+aJHSUKf2Pc8a/h1XLIwBsYMYmMGNrWMY1qI1jGIjWtaifZEa1ERE/SJ0EZ8gd0+R4h5fXF8yQ5+GzU+gi2mE5RlhcfOXlqSR6NBlp8kCrIprOK8tO+FBGCZEshmlEJKi/wmoZ16rljifF0lxrdNvstT01WFLOfYrdwyshhIVgnPayCSUdrjSTgYZQDKU3kr1Z4tcrSnIGCxnKEcOM5AzNNrM5MbMnOrbuBGsAAnwmMjxJ8RskZWxpscU6U0MkbUKP2qrXJ+4Nqvid7NoZBevKa14I1yLwhH2tqaC/8ArCPWL7or0URURjVERHtX2BIUTvwI5FCl+2nmu452tuVdfDgTB8VAvaap4svJ9MCqLp40CvNE1lnEIOVIPPp2X8dzKibLHFNJgEAd0WOr1E2rug9FQVOZrIlLRQxVtRXxwxK+titaKFXw4w2hjQ4UZiNFGixxNYIIBNaMY2NY1Ea1E6MdB//Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAOABUDASIAAhEBAxEB/8QAGAAAAwEBAAAAAAAAAAAAAAAAAAYIBAr/xAAiEAACAgMAAgIDAQAAAAAAAAADBAIFAQYHAAgRExIUFSL/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A6fuXdO9tt97d7Q6xzzadGttd472NNZPSegherGWNY2EV3KdSPoFUvsFjrrGukpVho1ymsWSRhWTk/wB4OQxgatid13rUyKJdj45suvDVLg99vmiOJ7nx5BRgZDJssXrRaHfHvrx9SrqifOmJQtjL4Fk6YyOjw6L6yPci7N0nrmlb88yLuO9it+jabeVlZ/IDCMHyhZ12yQSHdQfrJMGAiB52Vadd5uTwSmCpMdeirVyCjg4R5xGf5hgPM8DDGMvkWRRln4DOMf8AOcixHGI5lCPxCWY+Ai6NvusdLpsbFoG3LWtXI5QHYzUvAbiwMk4fW3XWwq2yr8xxDOF4NJjwwtgbIJEAQc5Hj8rUVSL1pZpVyKljdlVPcPLqgC3aGRUEgmWwYGOJnCKogCovNiZJBWEMA8xFCMcHgf/Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAvABUDASIAAhEBAxEB/8QAGgAAAgMBAQAAAAAAAAAAAAAAAAgFBgcJCv/EACYQAAEFAQABBAMAAwEAAAAAAAQBAgMFBgcIABESFBMVIQkXIjL/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A95GrOJBz9weI4SMwatNngmsnuhGDWMWaRCJXxxzO+ET2Ne9qN/8ACOVFT2/qW4Hz18YocbmWbnyC5xLs4aasg1z6wuwlDJvhhoR7N8T5AI3OhYf+RGOVqoqe3t7p/fTnahRUzeijIinIESjPQsYVXRvlE+gUhUMSR+7vyyDrKkcaJ/1KjE/iqipyR5YF4i6fEZ1cX/j66DuKgWi+kBeWnNeXLobwDOyiVTiDU0m8rroqxIn/AATkEmDRxTPWWRxLnuYkgdBOC9/yve7Tqv8Ar20r9JlcHqaeir9dVHyvD0TrbMVuiJeyOQSJ0K1RVhNUOb8pUe4NZEVvy+DT1kfhVjhsjW9bnpeSXXKKPQ7+Ozqq4wHKZmOwHjphBXMHy2L0OkoKZ9BLE/PHPHtZSLc6vnuDohzDJ4Iz0DpPPEsorKOqsx0nHJJAMIDmgK/X2ETHRvgKia9z2Ejvka9472KqKxEc32VfSuZ/BeStSJcWVlq+V6je0d8fDh9MZSkhrqufWJM8jMvvyRc/9qgIHa2rsmE5Bp32SqtsBU6wzvVyT908FOuTdd3fU8b17rC8921+Rrjeec32aZ/VU+jvUlnv7j43JoNGVmq+EZrB6wOae6mlIiaGK9UkRsJj/CnB7etzZdV5weVNeRck2AWazWt1dngdddLRTrBYsEyet/SaKwgCdEs/2m10zZR4VckskUi/MH08WMH0vnb+vp17Vc/vbfUdCdoqvN85hIHoMNXmVg7pKn9QtZWxVxdkUsl0QQ0RhV5Ka+7s3y2RxMjj1WPEbx91PjtYdwprbT6DaVWp3FBfZ3cbu5lutLphIcXT19hIbE0idayOqtICqgUZRw0mGDiJbFK2VJnnoL93Tyo4l4+QSN3+1rh72WeCMXDU0n7rZkRHMnkrZ0ztc4rQDwnTD/XZYTi/rY5JGNc1Fez0jG31vkh5x1JWSyHji/gOeIDiqL/qvXWGgdAyn3Sxb0Izm8gs+fNbEfJVsr7ciEc5HMNWFXxslcx2u0PgZU3XXes9x1W36Blt3puh6seoMzFzlrISHBMJDmoo2CaTJ6CMM2VWveXEz2SF0MKQNhRz/mwz/HjTuZCx3k75DMRGteMyK3561pETFa5HGe/OnI6T2RFT8aRNRff2an8RAkuA5aywOZKwt73DSdjvco2iqLUs6LKmmZkiDPgIlNLDnKUEsOGaP4nhs0Sl3BIk8JUpZDZfyvPUbwLjur5tsu8z3mp0OtB12yylzRabQWlMTpbUUXA0NYb+6SqpaoGF4NiMQBXsiAhc6uHGfK6aVXTPPQf/2Q==)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAcABUDASIAAhEBAxEB/8QAGQAAAgMBAAAAAAAAAAAAAAAAAAgFBgcK/8QAJxAAAQQCAQQCAQUAAAAAAAAAAwIEBQYBBxIACBEUExUhCSIjMVH/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A7Md+2TuhcumlV7bavWUzLJ0wl568bHbGFUBRT75kAja5kXPNhnW3wHJKpyWO+vxlinwb3vIavaqr3wM4GWeVTeep7BbxxinULW5XTziugeulEDn1HMzi+SeYYCE5VgjjDB9yzw/hx46md9UXebOaJs/Rm3B1mwYhTw4dYbGJl1quzyxyCWxPkI3TZdVn8JG4wWcQCcU+Gv48MAYByItaY79Wh4VDeQle2JLdJW7ZaxNJXmludXNxLZaZTjCVs8CwhuywdST+wrk6DwxzB6NK2zaEtW/S27VoqvbGiG8YKyGr8gp3VJk7hstY5Kvya2jcjgZfiIR2xW1R9UsoW2HLzlkuDqG0TqiV1o+2JZLZsqQ2psHYk2wf2ubwxHHREaGEG/bQcDB19MrL5hGEW0fnaoAqROoyRizlI/j45OgWvbetG/dL3KWvUOwbhLsdX6tpNKuDWo1Nx9FK2K22o84hrMv5VJ1lcjrY4IyWjXDAiQqliqy5FleEk38UTedHVxxIR89bNz1tjIRJHgLQ+LL7CrVeK2fqtL5pJmUVduWZ4iJJGwJ0woYlkB7lMg5znCc6RY9MUCQvVc2fmJWx2LFmj4pF2jTelZXtfAlwstXk5YQ/bfVqQLkRpKHKT1XhWrQhU8m6M40lOcJbvSISlGWxipGlCcJRnyvjnK0/0pXjz+7P5/Of96BGe1GMobTdnda+13Knf1i4SGotijcMZNb1iWUvVes85KuGXMIMt2rhyT5At1DSQSPCCpSpPjo6aGlahoWvrNsK21OEDEymy5dlK2sYOKWDuTiUPQjkAM0oSNu8efYOTSR8ZUt64Ukxc8k/k6D/2Q==)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAApABUDASIAAhEBAxEB/8QAGwAAAgIDAQAAAAAAAAAAAAAAAAgHCQQFBgr/xAArEAABBAIBBAECBgMAAAAAAAAEAgMFBgEHCAAREhMUFTIJFhcidLI0dbH/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A9ZW8OXPIHVO1idZQfFSY2MweLJS1TnK/e/YiercU6K07Mvocq3piTRlHi4mAlFFfByUElD5XuVltN+VnNjl+1rWt4Xx42Lx3Ld2FVw83Vqw4M+sOZYlVv0lKG4qPyymeU2l1DmXnU9o/PdpffHa7K412t2+vTVctQgZFXNjTIqZGOW2mPWI6sdWUvpdxhvyacbQpPktKWnPDOMrznPjXVyi4vcjds0Sv6wrdzq13p8Rsas3EObvBBtet9Sh6oiRHAjjpEQWxfqRJnCyryiJoxdaWhYaewjvy1ZYBqeNm4dv7bjbVI7U0FZdDqiz41mAjJuXRJFWBgtkxyQPwrEdHYawK8yOhSfF3yyTjPknx/cdMODIgTOSUsuRpKw3lIdSI+wdhr25Upv39lodHfcwjKnWHm0KQrGU4yrxzno6BE9gUUrkNyK2Hpq5bGtcBq6oaup0hihVOZXAv3oq+lTKpk20EjF4emR4D8rR+AW3g32xFyhKvNvL3Zckvxm2NGDMt1FMrubUrKcnTaLdcpST2tXBx1ISQLVHjQHWLa1IJeyYiMl5Wvjwo8YpkN4zJKsI5Tk1dtA6Ys9W2faKZPS+1lyUPBxRWpYvGdjmRiElqTFzUqE7HFIoin1N5kIaUkW42UI+I78MnIOVMYmw+Z7mt621Yb/xl3/DxJRTABEj9FqMmPiSfbeeFGbbi7idI9lNsEKZfWGhprCMoW40p1CVhgcWozXcfvDlnYdbGR8xXNgzWrL8s15TiMvytngrNKSj6MrY+QlD5JOVpFLZHIEx2adaQvKsYOpj44r0lcgbRuXUcAuLO227Bzl9hzwlREtH2AIUzCcTtfXhLcfPvKPNVMHs+9M0YhZeTDPDD2ToOK3FtS86Kv8jbxdDZv+tjKipVmuevnhCtljLgHU/SxrNGy78O1NQ+ESRyoRMVITUoO5iQwoAZL2FPQCn8UGmFODCscauTZjxL7bSGXKEB7ksO58MFvuEzbYiQ1u5ZStn35KSpba1jYS24puzEr/HC/kR39Cet8N9q/wDYPf8Ac9AruiZDb8vMbNv+zaLDavbt8zEiVOquixRl2bgqwzIgNH3ebryJQaQLPbJGejQ8zcmxDs/KFEUw04pKzpqQ/se/kv8A9ujoP//Z)

<a name="br2"></a> 

2

Fast Approximations and Coresets for (k, ℓ)-Median under Dynamic Time Warping

properties towards the problem at hand and can be used as a proxy to run an algorithm on,

producing a solution with a relative error of (1 ± ε).

Clustering and especially k-median represent fundamental tasks in classiﬁcation problems,

where they have been extensively studied for various spaces. With the growing availability

of e.g. geospatial tracking data, clustering problems for time series or curves have received

growing attention both from a theoretical and applied perspective. In practice, time series

classiﬁcation largely relies on the dynamic time warping (DTW) distance and is widely used

in the area of data mining. Simple nearest neighbor classiﬁers under DTW are considered

hard to beat [[31](#br25), [44](#br26)] and much eﬀort has been put into making classiﬁcation using DTW

computationally eﬃcient [[30](#br25), [36](#br26), [37](#br26), [40](#br26)]. In contrast to its cousin the Fréchet distance, DTW

is less sensitive to outliers, but its algorithmic properties are also less well understood, owing

to the fact that it is not a metric. In particular, the wealth of research surrounding k-median

clustering for metric spaces does not directly apply to clustering problems under DTW.

For time series and curves, k-median takes the shape of the (k, ℓ)-median problem, where

the sought-for center curves are restricted to have a complexity (number of vertices) of at

most ℓ, with a two-fold motivation. First, the otherwise NP-hard problem becomes tractable,

and second, it suppresses overﬁtting.

The construction of ε-coresets for the (k, ℓ)-median problem for DTW is precisely what

this paper will address. To this end, we adapt the framework of sensitivity sampling by

Feldman and Landberg [[25](#br25)] to our setting, derive bounds on the VC dimension of approximate

range spaces of balls under DTW, develop fast approximation algorithms solving (k, ℓ)-median

clustering, and use coresets to improve existing (k, ℓ)-median algorithms, for curves under

DTW. We rely on approximations of nearly all objects involved in our inquiry, thereby

improving the bounds we obtain for the VC dimension of the range spaces in question and

broadening the scope of our approach.

Our analysis of the VC dimension is possibly of independent interest. The VC dimension

exhibits a near-linear dependency on the complexity of the sequences used as centers of

the ranges, yet it depends only logarithmically on the size of the curves within the ground

set. This distinction holds signiﬁcant implications in the analysis of real datasets, where

queries may involve simple, short sequences, but the dataset itself may consist of complex,

lengthy sequences. Note that our results hold for range spaces that are deﬁned by small

perturbations of DTW distances. This means that for any given set of input sequences

requiring DTW-based analysis, there is slight perturbation of DTW with associated range

space of bounded VC dimension. This is suﬃcient to enable a broad array of algorithmic

techniques that leverage the VC dimension, particularly in scenarios where approximate

computations are allowed.

Related Work Among diﬀerent practical approaches for solving the k-median problem, a

very inﬂuential heuristic is the DTW Barycentric Average (DBA) method [[38](#br26)]. While it

has seen much success in practice [[1](#br24), [28](#br25), [39](#br26)], it does not have any theoretical guarantees

and indeed may converge to a local conﬁguration that is arbitrarily far from the optimum.

Recently, theoretical results for average series problems under DTW have been obtained.

The problem is NP-hard for rational-valued time series and W[1]-hard in the number of

input time series [[12](#br24), [20](#br25)]. Furthermore, it can not be solved in time O(f(n)) · m <sup>( )</sup> for any

o n

computable function f unless the Exponential Time Hypothesis (ETH) fails. There is an

exponential time exact algorithm for rational-valued time series [[10](#br24)] and polynomial time

exact algorithms for binary time series [[10](#br24), [42](#br26)]. There is an exact algorithm for the related

problem of ﬁnding a single mean curve of given complexity for time series over Q<sub>d</sub>, minimizing



<a name="br3"></a> 

Conradi, Kolbe, Psarros and Rohde

3

the sum of squares of DTW distances to input curves, which runs in polynomial time if

the number of points of the average series is constant [[18](#br25)]. Furthermore, approximation

algorithms were recently developed [[18](#br25)], and some of these can be slightly modiﬁed to work

within the median clustering approximation framework of [[15](#br24), [17](#br25)]. Unfortunately, known

median clustering approximation algorithms either have running time exponential in the

length of the average series, or a very large approximation factor.

Approximation Algorithms for Series Clustering In the last decade, the problems of (k, ℓ)-

median and (k, ℓ)-center clustering for time series in R under the Fréchet distance have

d

gained signiﬁcant attention. The problem is NP-hard [[11](#br24), [13](#br24), [24](#br25)], even if k = 1 and d = 1

(in these works, time series are real-valued sequences), and the (k, ℓ)-center problem is even

NP-hard to approximate within a factor of (2.25 − ε) for d ≥ 2 [\[11\]](#br24)[ ](#br24)((1.5 − ε), if d = 1).

For the (k, ℓ)-median problem, all presently known (1 + ε)-approximation algorithms

are based on an approximation scheme [\[16](#br25), [22](#br25), [24](#br25)] which has been generalized several

times [[2](#br24), [17](#br25), [32](#br25)]. The most recent version of this scheme [[17](#br25), Theorem 7.2] can be utilized to

approximate any k-median type problem in an arbitrary space X with a distance function.

All that it needs is a plugin-algorithm that, when given a set T of elements from some

(problem-speciﬁc) subset Y ⊆ X, returns a set of candidates C that contains, for any set

<sup>′</sup> ⊆

T

with roughly |T | ≥ |T|/k, with a previously ﬁxed probability, an approximate

′

T

median. The resulting approximation quality and running time depend on the approximation

factor of the plugin and |C|, respectively, with a factor of O(|C| ) in the running time.

k

For the Fréchet distance, plugin-algorithms exist that yield (1+ε)-approximations [[16](#br25), [22](#br25)].

For DTW however, the best plugin-algorithm [[18](#br25)] has runing time exponential in k—roughly

with a dependency of Oe((32k<sup>2</sup>ε<sub>−</sub><sup>1</sup>)<sub>k</sub><sup>+2</sup>n)—and approximation guarantee of (8 + ε)(mℓ)<sup>1</sup><sub>/p</sub>

with constant success probability. Here, the Oe notation hides polylogarithmic factors. In

principle, some of the ideas from plugins for the Fréchet distance could be adapted, but the

more involved plugins, i.e., the ones yielding (1 + ε)-approximations, crucially make use of

the metric properties of the distance function.

In practice, an adaption of Gonzalez algorithm for (k, ℓ)-center clustering under the

Fréchet distance performs well [[14](#br24)]. Similar ideas have also been used for clustering under (a

continuous variant of) DTW [[6](#br24)], but there are no approximation guarantees, and the usual

analysis is based on repetitive use of the triangle inequality. To the best of our knowledge,

all (k, ℓ)-median (1 + ε)-approximation algorithms for Fréchet and DTW are impractical due

to large constants and an exponential dependency on ℓ in the running time.

For the Fréchet distance, ε-coresets can be constructed [[7](#br24), [19](#br25)] that help facilitate the

practicability of available algorithms. Using ε-coresets, a (5 + ε)-approximation algorithm

for the 1-median problem was recently analyzed [[19](#br25)], yielding a total running time of

roughly nm<sub>O</sub><sup>(1)</sup> + (m/ε)<sub>O</sub><sup>(</sup><sub>ℓ</sub><sup>)</sup>, in contrast to a running time of n(m/ε)<sub>O</sub><sup>(</sup><sub>ℓ</sub><sup>)</sup> without the use of

coresets [\[17\].](#br25)

For DTW, no coreset construction is known to this point. This is at least partially due to

prominent coreset frameworks assuming a normed or at least a metric space [\[25](#br25), [33](#br26)]. Also,

recently a coreset construction relying solely on uniform sampling was developed that greatly

simpliﬁes existing coreset constructions [[8](#br24)], including the aforementioned coresets under the

Fréchet distance. Unfortunately, the construction again relies on diﬀerent incarnations of the

triangle inequality, limiting its use for DTW.

Results To construct ε-coresets, we use approximations of the range space deﬁned by balls

under p-DTW and bound their VC dimension. Assuming that the input is a set of n curves

![ref2]

<a name="br4"></a> 

4

Fast Approximations and Coresets for (k, ℓ)-Median under Dynamic Time Warping

Figure 1 Example of a traversal between the red and blue curve realizing the dynamic time

warping distance. The sum of the black distances is minimized.

of complexity at most m, we present an approximation algorithm (Theorem [40)](#br22)[ ](#br22)for k-median

with running time in O(n) (hiding other factors), that improves upon existing work in

terms of running time, with comparable approximation guarantees. Our approach relies on

curve simpliﬁcations and approximating p-DTW by a path metric. This allows us to apply

state-of-the-art k-median techniques in this nearby path metric space, circumventing the use

of heavy k-median machinery in non-metric spaces which would incur exponential dependence

on k and the success probability. Our main ingredient is a new insight into the notion of

relaxed triangle inequalities for p-DTW (Lemma [22).](#br15)[ ](#br15)We then construct a coreset based on

the approximation algorithm. For this, we bound the so-called sensitivity of the elements of

the given data set, as well as their sum. The sensitivities are a measure of the data elements’

importance and determine the sample probabilities in the coreset construction. We construct

an ε-coreset for (k, ℓ)-median clustering of size quadratic in 1/ε and k, logarithmic in n, and

depending on (mℓ)<sup>1</sup><sub>/p</sub> and ℓ (Corollary [42).](#br23)[ ](#br23)We achieve this by upper bounding the VC

dimension of the approximate range space with logarithmic dependence on m (Theorem [13).](#br8)

2

Preliminaries

ꢀ

ꢁ

We think of a sequence (p , . . . , p ) ∈ R <sup>m</sup> of points in R as a (polygonal) curve, with

d

d

1

m

complexity m. We denote by X the space of curves in R with complexity exactly m and

by X the space of curves with complexity at most m.

d

d

=m

d

m

▶ Deﬁnition 1 (p-Dynamic Time Warping). For given m, ℓ > 0 we deﬁne the space T of

m,ℓ

(m, ℓ)-traversals as the set of sequences ((a , b ), (a , b ), . . . , (a , b )), such that

1

1

2

2

l

l

a<sub>1</sub> = 1 and b<sub>1</sub> = 1; and a = m and b = ℓ,

l

l

for all i ∈ [l −1] := {1, . . . , l −1} it holds that (a , b )−(a , b ) ∈ {(1, 0), (0, 1), (1, 1)}.

i+1 i+1

i

i

For p ∈ [1, ∞) and two curves σ = (σ , . . . , σ ) ∈ X , τ = (τ , . . . , τ ) ∈ X the

d

d

1

m

\=

m

1

ℓ

\=

ℓ

(p-)Dynamic Time Warping distance (p-DTW) is deﬁned as





1/p

X

dtw (σ, τ) = min 

∥σ − τ ∥<sup>p</sup>

.

p

i

j

2

T ∈T

m,ℓ

(

)

i,j ∈T

ꢂ

ꢃ

1/p

P

We say

∥σ − τ ∥<sup>p</sup> is the induced cost of T.

The central focus of the paper is the following clustering problem.

(i,j)∈T

i

j

2

▶ Deﬁnition 2 (Problem deﬁnition). The (k, ℓ)-median problem for X<sup>d</sup> and k ∈ N is the

m

following: Given a set of n ∈ N input curves T = {τ , . . . , τ } ⊂ X , identify k center curves

d

m

1

n

P

C = {c1, . . . , c } ⊂ X<sup>d</sup> that minimize cost(T, C) =

min dtw(τ, c).

k

ℓ

τ∈T

c∈C

![](data:image/jpeg;base64,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)![ref3]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAWABkDASIAAhEBAxEB/8QAGwAAAQQDAAAAAAAAAAAAAAAAAAYHCAkCBQr/xAAzEAABAgUDAQYBDQAAAAAAAAABAhEDBAUGIQAHEjEIExZBUWHRFCMyM0JSU3GRoaKx4f/EABkBAAIDAQAAAAAAAAAAAAAAAAYIBQcJCv/EACURAAEEAQMEAgMAAAAAAAAAAAIBAwQFBgcRIQASMUEIIhMUUf/aAAwDAQACEQMRAD8A6itzdzr/AKbuFdUlJ3fXqbTKdWIiZWTQsKSpKjFBQg96kBA4hsDqceWkTC3Z3JCAvxvcEMxXi8QtvrC+fnnyXf0/LOsd2wIe5l5fa7urz3ErdbKhxEBH0n6c1OSQ+D5Bm8KeKlhwwJ92z0f0YdPIAa5YtT9TtQo2oOXomYZWIjlk5BFLEdkROOPttwi+NuOeOUTp/MVxXGX8Zr3nq+M465Gjmbhx2iNTJoCLciBfaqvKr59+3coW5+4M7cFvy0xe1XjS0eu0eFMwZlPeJjQIlRlkqhJaOSgqcAqbCeQbJ1bPzh/iK/l8NUrWsCbotsl3FwUUukkAkVKVGfVurH39H1dE6/f9P81oL8B87vbXEM7k2k+daSUvqgP2bF788hQSsJRDvRV+o+BT1x/V6pXXKrra2wx1mtjtRmSr5ZkLLbbXeavtIpH2CiEW2ybrvsnG/HUML47K9fua7a9cEvddHloFXnJqahy8eSnVxoKZhaFJStaCEKUnhnjglmZtJY9ju5ipR8ZUPJJzIT5OfXIGQA+G/vRo0X5l8cNFLLIrqbNwWG/Jk2smS+6trkAK4+a/ZxRbtgBFXfwIiP8AB6jKvUHL4dZGixrgm47bLQg3+lWmgiIAgohHDI+ERPJKq7c9bKidkW46dWaTUIl30SKiQqchOqQmRnwuIiUmoMwpCSVcQpYhlIJDAlyWfU+vkJ++P3+GjRphPjvpRp9hNJfw8WxxipjTLGG/Jabm2khHXW434wNSmTpJiogqpsBCK+VRV56B82yS7vpUB22nFLOPGJplVYjM9gE42RCgx2WRXdfZIq+kXbjr/9k=)![ref4]![ref4]

<a name="br5"></a> 

Conradi, Kolbe, Psarros and Rohde

5

X<sub>1</sub>

X<sub>3</sub>

X<sub>2</sub>

Figure 2 Illustration of a coreset (red), i.e. a weighted sparse representation of the original set of

curves (in red and black). The weights in this case are w(X ) = 3, w(X ) = 2 and w(X ) = 1.

1

2

3

An inﬂuential approach to solving k-median problems is to construct a point set that acts as

proxy on which to run computationally more expensive algorithms that yield solutions with

approximation guarantees. The condensed input set is known as a coreset.

▶ Deﬁnition 3 (ε-coreset). Let T ⊂ X<sup>d</sup> be a ﬁnite set and ε ∈ (0, 1). Then a weighted

m

multiset S ⊂ X with weight function w : S → R is a weighted ε-coreset for (k, ℓ)-median

d

m

\>

0

clustering of T under dtw if for all C ⊂ X with |C| = k

d

ℓ

p

X

(1 − ε) cost(T, C) ≤ w(s) min dtw (s, c) ≤ (1 + ε) cost(T, C).

p

c∈C

s∈S

▶ Deﬁnition 4 ((α, β)-approximation). Let a set of n ∈ N input curves T = {τ , . . . , τ } ⊂ X<sup>d</sup>

1

n

m

ˆ

ˆ

be given. A set C ⊂ X is called an (α, β)-approximation of (k, ℓ)-median, if |C| ≤ βk and

d

P

ℓ

P

min dtw (τ, c) ≤ α

min dtw (τ, c) for any C ⊂ X of size k.

Relaxing the problem to (α, β)-approximations allows us to pass through so called

d

ˆ

p

c∈C

p

τ∈T

c∈C

τ∈T

ℓ

simpliﬁcations of the input curves.

▶ Deﬁnition 5 ((1 + ε)-approximate ℓ-simpliﬁcations). Let σ ∈ X<sup>d</sup> , ℓ ∈ N and ε > 0. We cal l

m

an (1 + ε)-approximate ℓ-simpliﬁcation if

∗ ∈ X<sup>d</sup>

σ

ℓ

inf dtw (σ , σ) ≤ dtw (σ , σ) ≤ (1 + ε) inf dtw (σ , σ).

∗

p

ℓ

p

p

ℓ

σ ∈X<sup>d</sup>

σ ∈X<sup>d</sup>

ℓ

ℓ

ℓ

ℓ

A range space is deﬁned as a pair of sets (X, R), where X is the ground set and R is

the range set which is a subset of the power set P(X) = {X |X ⊂ X}. Let (X, R) be a

′

′

range space. For Y ⊆ X, we denote: R = {R ∩ Y | R ∈ R}. If R = P(Y ), then Y

|Y

|Y

is shattered by R. A key property of range spaces is the so called Vapnik-Chernovenkis

dimension [[41](#br26), [43](#br26), [45](#br26)] (VC dimension) which for a range space (X, R) is the maximum

cardinality of a shattered subset of X.

We are interested in range spaces deﬁned by balls by the p-DTW distance: We deﬁne

the (p-)DTW ball, of given complexity m ∈ N and radius r ≥ 0, of a curve σ ∈ X as

d

ℓ

B<sup>p</sup> (σ) = {τ ∈ X<sup>d</sup> | dtw<sub>p</sub>(σ, τ) ≤ r}. Deﬁne the range set of p-DTW balls as R

p

\=

r,m

m

ꢂ

<sup>m,ℓ</sup> ꢃ

{B<sup>p</sup> (σ) | σ ∈ X<sup>d</sup>, r > 0}. The p-DTW range space is the range space X = X<sup>d</sup> , R

p

m,ℓ

.

r,m

ℓ

m,ℓ

m

3

VC Dimension of DTW

In this section, we derive bounds on the VC dimension of a range space that approximates

the DTW range space. Our reasoning exclusively relies on establishing the prerequisites of

Theorem [7](#br6)[ ](#br6)below.

![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCACFAhIDASIAAhEBAxEB/8QAHgABAAEEAwEBAAAAAAAAAAAAAAkEBQcIAQMGAgr/xABLEAAABgEDAgQEAgcFBgQDCQABAgMEBQYHAAgREiEJExQxFSJBURZhIzJxgZGh8BcYQrHRChkkUsHhJTNi8SaS0ihDV3KClcLV4v/EABwBAQACAwEBAQAAAAAAAAAAAAABAgMEBQYHCP/EAEURAAIBAwMCBQAHBAcGBQUAAAECAwAEEQUSIQYxEyJBUWEUFSMyQnGBB5Gh8CQzUlNiseEWJUNjwdEXNDVE8UVyc4Ki/9oADAMBAAIRAxEAPwD9/GmmmlKaaaaUpppppSmmmmlKaaaaUpppppSmmmmlKaaaaUpppppSmmmmlKaaaaUpppppSmmmmlKaaaaUpppppSmmupYxylA5TAAFMBjgJeoTkAB5IXuHAiPA9XfgAHt31bxer8BymYDcmV4SKC3KJRABSHkSdCo9ZRL2NyBT/bSlXXTVtTeHOJgMQUxBVVMhVAAgKFTN0gsQ3cDJqcgZMewGAOedfR3SqQGKJQUP+qQR/RkMYAETcjwYC8h3J+t1AA+3AjpSrhpqhZuVHAAcSm8swdQCoTyjlA36gdHJuoBDkRMIhxx7Dz2rtKU0000pTTTTSlNNNNKU0000pTTTTSlNNdSpjlAOjsIjxz09QF+vJu4cB24559xDsOqFR6cp0gKQVCqKlSOKYAYqACQ5/NVOPHSmAk6BHgfmOUPceNKVc9NaO+ITk/c1izZ3nm87O6V/aJudrVUbyOI6KlWlbn+IJk1ihWSzY9fRexi0ib4M6kXvpiOkhIKHPWYCDzZtge620bjcGVNLM1el8bbpqTTqQx3GYltlYGkW2p5Al6+k7eyoUtaQk14yqWV43kJWonGQc+rgxRWUOkbgplK3701a2LpwscSL9JOkolIIl6TuugQKd0VPv5SJjcCmXqP1FOA9QccDdNKU0000pTTTTSlNNNNKVTOgJ0k5APM6/wBCYwc9CnSb5g79hAvVwP01ZjmOBSAAqpm81QiB3CAqnB2Am5XTDzCfoTpFVAg9QdQGKbn6DfVUirF6TGMAd+ekeOQEolEB7DyHA88fcAH6axplbGDbKOOLDjs92yHjxCdj2zH8ZYus34Sv0Em0dNXYL16x+hkPhjlb0gN3C3o1vNaLuUekPN6iqV69IxjFJ5fndIrCBTIpeUmmuUD+oVch1G6SKqgcSjz2E5S9+edDgQPO8lMwcHAy5Wxf0qqom45KXkBHgxhObgexQMP2DUVvkb/tpip7Tbb833v7b4gx4M+OqLiVzX90tQpyTgharbUJwbnZUtxF4bA3iYKzRnwrG4TRJeYvQPmPwcIB98vPGS2ErybaowOYvxlnsjgqBNpdTYoS+6EbgigdGdxwni9WRZH/ABxVR+JHscJ8X4jQg5cCuXAtQFSQiyZDY4BIzgc5HbPcn29akQyykNG5TwT4zkHGUT7wPpj8+/bPNZY8QbxGdv8A4aeLK7mHcsnf1ahbru2x6y/s4pZbdJlm3sHM2BqrKNzy8URnGEi4J8ku+MqcoOBQTBPlXtCXsC8ZTev4r+9+sM9qG3es438PHEM/b2OfMmZGVdTFivsQASSOPU6VIJxsW0qtxkTBES7+lGNOgWObzCwTQjHAVeU6d2iWnxB4qbtm9RfIVIw7a2zyLp2zyIsYR9WXx68aLpR0puJjXMWsWz5HkXfwe9IV1JtE/wBk9laFrHxS3Fi/jDrM+Pqts98LvbjR8N01Cv4hxZCdNWqEU5BAJy7XBdBdy5eLFQSQNPXa6v27uUevDEbhLzD1ZcAa+f0BJaOEhSVTB3EdgcYOGJ7Z58pPOfQnNTZWZ1W9EunRPHM7CHxgm4S7sKxslBBhdTkl/MDg4HAatyrTcoesVidtEu5WUiq3EzE3KlQAFnZW0Q1cSD1JqQTkFVwmi2UTSbiJAE4AXzA9w/Nbkjxgt2W8bJUXhLwpcaQrexJi2schkrPlYcOoROohCLpv2spX2btIK/02d5DoNZwZl91pF8v4eAuxMjK2xx5nPdSqsbMadi294qAhVoKpY5uCiGS77HSJPNYqZSdKQzZOrxIwSi0XZcYg1mvMlHaSoWcAjPLd5e2nbCdvmztTIq2FoKUhEco2ALJOR0g/ResY9ymR0iiwgUCM2xouJTRcnInH+YuBAIiHnCKYCPjr+HqTVNQin068W30lZHEsQBSSVQowI9wXZubd4hK7kXBQmvufRt9+x/o3pjra0680ibr79okdtbR9DWFvdqek9KZpT9Lv+ontXMmoX9tA0TWttBKsHjLNHctz5cm4qkc8FxfjcuXFaWOVy0KnlycNcWVCvDkIK9HBdBgQOAnCGGyfEhiwMImBj5AGER501msscwSKVNMrghEygQhCGECEIQOkpShx2KUAACh9AANNenFuwABjGQAP6kHkbfXdz6c4/wCufz03VMpZj9VaYMsTgXUgAyRgAeHwBwAPQflx7PTTTWxXQpppppSmmmmlKaaaaUpppppSmmmmlKaa45D7h/HTqL9w/iGlK51wI8AI/bTqL/zB/ENfJhDpHuH39w+g/wDbSn64+fb5/SvgFgEeOB199Yfn/X79UfUA/wCMeBHj8v2frDz/AF213gYvHcQ9u3tzz+ft/D+WpkViFKeXOM9z7Y9z7g1SMoOPHWZvZcA544/n19jX2CnI/qmAOO/Id/3d+NfXWH2H+X+uqQDlHqEohz+sHtxx9eB5HnkB+mnmD/6e/t+fHv8AXVcOvBDMcA5VTjnt35z/APNGbaQG2Rnjyu+DjjngEc5J7/8AaqvrD7D/AC/110qO0kjlIfrL1mKQhhAAKc5+ekhR57nHpHt+WqRw9QaoqOHS7ds3RKJ1l3CpEUUiB7mUVUOUiZQ5DkxjAAc++vMOLTVlFDl/EkGkIlKVRQkvHkU5DkU1E1fVD0ij8wG+URDrD97J/sOe2cD/AF+fk/FR4i+skI/Nz8fHrk4/Tn1r2IO0RMcvJgEhwIPUHHJjc8AXv354H+WgPEfkARMUxwEwEMAAcADjqExeewF5ADfYRAPcdRs3/wAVzw6sb3G247v+8HCVZuFUmX9csFdmLGJ5WInY1YUZOLctk25iCZBUCk+RcQAQMACIDryhfFf2sXIqKm2uKzRvNbRpRG0q7R8dmy0hj9056Rho66iMzAfApGfSSkVa6kHq/WIxEyImJ6Xg4tyR4b+mBgZ+e/8Ap8Zqdy8fawc/8w8du/GP7XYn0qU4X7coAJhOXqOJCAJeBPwAj1l790x47GHjnkPvp69uApEOJ0lFi9RE1C9JxAA+bkORAOgRAD9/lEwB9dRXOPERyNa0nUHiLw+t7gZGmUfR1QudMQmwtihWVKIC1NfcsElrqekQooA4A8sWpzgkWMij6Q3m9ZbAbcF4uQHIP+7t23AUQFEiTrfs+IVUXAlOQoJf3ahA66QJDyTqDyeoxRE3UAgBz+CT9APj9ex+c/lTevH2sHOM/adjkA9xz69vjGal2K4SOUpim6imADAIdwEDByUff6hyIfcAHXWD5qYDCVUDdJwSEA9/MEDCBADtyI9Ju35aiRSxJ4ulzVXtyW7zbZgotiUPKqYkNtaeZkPj0jw3mp0suU/7ZaF+Pk60Bxjm9v8AwbWRnCFM++BR3X6Ytchsg3YZPQXU3U+INmF2+iRK3piOzGNV2oMDRrnkZkchsFbFlkcgSorIRZoF6VeA+Btwl0RRffE+tsBP93J+4fHH+fv6VO5f763P5SH4P9n8/X2xn1lWUmo5BJRdw4TboI8g4XcHTRRbHL2FNwqocpEjh37GNx8pu/bVnG900BJzZoIAUOCaYjLxwAc5gExClH1PcVClMYn3ApvbjUZgeFxQpkDw2T90G9/P2OnPKFkw7m/cKN1xRfmKXY0Feq0FPjRlYgynQuePB6iKiqKR/OL5fBrJO+EZ4TtaRO5n9m22uIRQbOXi7+Vhk2xV0GqQrOVTpLP0x6RSIZZU5DKdAJ8cfMHLJ/u5f3D4z7/OOf8AWAy55mt+wP8AWfln09M9gCR84IGYcg+LV4beNLhYse3veRhKr3anTz+u2auyloFGQhZyJWVayUa9Ik1WT9QxckMi4TIoYoHDscQDvqHue8cPbTjDbtlDNW2OrW/fEGL6yvZrm92+xas1jGjShnDdCPY5evxSLKUJKyprSLiDXJXp0X6UTIm8lL03zewofiHeGtA02BxNtBPFbk1aRBRUPWsAbPKI1ybfazQ4Nr6FabaVhdxEH/CEAc8THST48gJmLiTi0DEWF15id3JkjxH8vsnVfxXs9wvtzxhkFMv4SzJlvIITORKXXFf00PZMk7TUaDBILWdJmYEZbGh8vkLHvHSqYWhf0XLgGz2jkz29PcDPr78+wyfapLL/AHsPpj7Q8j48v6D3P6458IHxTqf4n+3ORysFQg8K3+vXWSqtowepaPxLP1OFFRwekysg9dw1YWM2skdHyrliZKKEDJsVjCYeNbe7g9+Gz3a69q7LP+4fGOKlrm0lHlUJbJ5FMkwhCrMUZFdAzUrgOWR5FmmJD9IiDnkoiAajitHg0xm4W8Y4zfu13M5Ys25fE686njS6bdCl284eprIUkWbNu1xOm/vhkJiQQ6CWpY9uV/EaKBykJGdAmN5nbzihbwm0bS1y9gDbU224WU8fI5H3jYGqBsaMIqYEjo1aPlLAD1a7BTcfUpmvNxr7KKeT5XzpCXikBqyPxDqb03SLy8THgHC8nuuc8Aevv/DtDMTzslgXI8zRlxj8/L3/AJB9KXZb/tDG1rdTvdzbs+kIcKCtD3uXrO3LMCEypP483JwUBIPGJZKPmXEPAo1qckmxo+SrlfA84SYj1pZx8Ua/CQI7/QelJsVkhXScJqtwTKqDhM5FEDpKceUqmqmYxFE1gHqSOQRKcvcB41q+jV9ru7vCgkbQuN814Py3UY9dQGzWOl63dqjNpov49w4K3SbuFmT5JNFy3Axmp1S8mAC9IhrURXD28TZq6K22i/gHOe2eKOMZWtpl1Wc43sWJ4lQQdPP7Ocyk/GZZGnVNlHhB46wuOPmApNX5ERupfRf8TZW3jdtKewPfsOaIVK+WTxccFinhnI7grk9vf+AqV/1iXPYqgh0ibqAodPYQAS89X6wCPcPyEeeA1UkOBw6i88ciAc/X8w9+w/QdaIYN34YHz1dFcPoLW/FGeW8a+tMxt5zxWgoWcoWptnTNqys0jSDyMmDeuSgyLNWHkySap3KC6SgNSgceneZnyCBQMbk/UbqATdRim57lMYeOsxR5AxuC9Q9+kB5DVqvVXprjkPuH8Q05D7h/ENKVzppppSmmmmlKapxXAoiAh3AOR4Afb8h/r/LVRqwv1AK6KQAVMoKCpkik/wDKMQDEBTzePp3Dj8vbvqRjnccADOf1H+fb9ffFVLeZE7GVgiuRlY2YgBnH9nJx3HfvVQ4kWwFKU5ukT+YPQYekRIQQKJ/ryTrMUA+/UH7NRSeLF4jy3hrYZxVlscRkzGwvucaviSagizq1YGDhZmv2iwvZwj1CDsRnblmaspIs447RsV66cIt/XNxOB9XHe94q+2Tw/wDIuHsZ5xjcqu5zNjCUeUlPHdLJbI1q2hZ2u1x26lTmlow0emSQs0aQ5yJuRKmdQ/QIEHUlrmLjnqZEHka0kE+sVARcoIO0+sCmKVcCuEzFTEpTCAKgQTAY5Sf4xEMM0oAZYiC4HCg7iCcbcjPbJ/8A25weK2rS3ktwJdVhYRSSNFbICYjcm2YLcskhBypLKEZVcA5yM5A1v227x8C7sKg3uOIrYEk2Mzg3sjXnQpsLnWZCeYC/TrtohFDnPFzrAnnoPowi7v0yzZcgrD5Q9VBuK2mVPPsvAXqHuuSMO5jpkXIRFay/iCxlp9wPBOV279zRrTKEZPDWLHkrNxsNMWKpHKzNMuINgBJNkVIRNectbVcO5rXh5ifgFq3kCpldjRMh1B0lC3ukqLKFI5lavMotDBGSjlATtHSp2zgSoOnKPA+Z1BgFSM3zYJV+GV1lTd0dEScgZi/s9lcY0yHVKvHAdNBpILhGXUuWbc8jhKo5mAJTyvJFFRP0SQP+pvrwtdKN128UfbjYwB55xgu2c8YAcdiSucD0x0fRtXAk0bU00ycgbdL1uRYfFI2jbBqjLHbSMwyw8aO1YndGiyMEMl7w/u7tdJutW29bv6XP0LLL94eo1fNzCumjtuGeJxNuu9rydFsq7wwweRrjWo2St7/ExkpYtKRj5qELbrAMSEg8kK+Ks+tRIVSkUTIVQxT9h8s4gAGD8uTFD9oh9O+o/kcwbbN3cJcNu+TYeRgLFLQnrLThzKMOlTr1BwrWSj1q5ajw6zmQBg4RfDEylakUnR1BU9M5BISgKY4nTnsz7BHyrC2k3D7xNtEyHqIu9KnPmLcXh+xEL5sv/aCp0wadsxVINBl5wboU8COOW8ZGVj4JZRmRlmOwLq2LbDuDkgckY9O5zjJ5GMg8HsQRXA1DSdT02cR39pc2shQNFbCPxDcI20pPHOdq+BIGDRvjDKQwyGAMrhZBIRL1dSYGOJQ6y9PVwI/MX5h5KIByUfqHHbntqoKt1gIlKYoAYxfnKJRHpEQ6i9+5De5R+oDzrFVBvtMyPSK5kWgWWItNKtUXHzsBOwhyPY6cgJBumtFO2K5VC9JF2zluqUTF6yEHg5OoR4ya1MBkvm45KcxOCCAgHSYQ4HsHze3Idx/PWyYyMNnKnGBjBOcep+COcY5rn/bjLTW30ZQcbHl3y5Jx90LjA53HJx81VdY/UA/r9+vvqDpE30D7/wBfnrqDp+vPP5CH7tcm7cl57fnx9g/LVGGPQ49vU++D2qodXyEODx3Gcds59+9cAuUf8I/v/wDbv+7XYBwEOQ/Pjn8gEfv9OO/7Q1RFEoiIAB/48/w5D/T+WuwTmAB7lLzz78CIfl9B/d9PrrGZC2AkLgnH3jxzj/qT88dj6UjJ7iVpx7JFtPp8nPx6nnjmuwy/TzyHsAj2AfoHP3/cHPHfXWLxISc8HARAvJTAACQTk6iAp3Hp55AvPfgwgHfVhsU1EwEHKz1hkGsbBQzB3JS0i8VIgyYsGbdRd06eKHMBSt0USHOfkQEeA6R541F+7ytlXfcaOquBzZh28bYiplcXTccEYtRMiZgh1UBc1BntWeA4fEb1ayNxYW0+ZVFXnkwrA1SClKjbPjURbzKu+QADPYZJ/j/Prz2q5LSYWJH8bv4TYBZePXPBPccEEZ/KssZt3dy0PcZrbvt1xpdMzZ7H0zV9N12IAcM4TmbCVOVjHGdb6iq8GoKGraj2xwMQWHfmthmTaHO6hhkjPWkBPiAf7MLIbpMssN2dR3xZAq258K4Sfv1xs9cCfirrlGDQava2vQUELNCLYkosY9ZljoKvHd3JSuwwskTykl8PMZ1+hJu620bDsUVyqGn2VJjJ2ceC2fTr0jy85dyZYiP52xys28KkirN5Kuk18Um5Z6dJMstNunqgEa+pAhPDlqmfdyogtmpWxbfsSMwM5rVRx/c1EsiZHj34CuwNk6SUh2SdRYhAKOI6142BnNgs/cC5/EyQRflOtZ7qOUBLYnxgwLow83hjG4pk4P6kcHChjXfs+nnuoYr29leytA+WUtsllKAF1tIuJL1VbyyhRtiJDSMq7mGsW03c9uutmAsd4fTrNDytuoxzFsMdbiLTE5lc2fG9bsFSbFqr64yt7NR2Skxka2v2g3J9joYlmVYFJkPxSYI8FnOxFrk9r2xGOb7g95+4qJip+0LkpCGT83TqbaBi5OVTPZJqo48auE3JKjXZWRhVZtnAGdyJmiMa2R+ILCh1qV1y3jbfNv8APBtywrUJrNObKzA11MMB7fIFG1Wah1p5DNi1WzXkWrhAlbxoidxBxMhaxCSUik5ViYY1yY4p61hzfsT3N+JJiK9403w3mI29Yav0O9I3207cJs9kk4a2xrhNOKsl5z+6aQhcnU6WfJfi15jwmMqmpGyvwxqFndhEi4e3W1+kqFY7cHB5IDFcEbvg/wCIk5HAFbN9rscobSOnoYbK3Chbu4mXxJrlQBuWacBdhBBCxQoqgMVkeQ4atQt0X+097NsY5ZoWGNrcFLbvbjaMhY9qFjt1LlfhmH4GqXxFIC2WMvzOLsX4gnol27jY57W1IeJAizl2PxUPRdK/6eI9E5UETAUxSqpEMJVBETp9RQN0HHgOTk54MPHvzyAewRq+Gp4Xu27w0dvULhrD1ZZTM1IDHymVsn2SNZjZ8jXNo3ORxZnQHSXUaxRHSzolZhPVOTQkOu2YfEXotvUKSbsjdSHPSYoeYoAdSXk9QAcQ6yl6j8kP+sQ3PzlEDcBzxrJtIYFHaNAADEhCxHByMoBj37Y9686rCPxRGiRrK291UEKX483OWz37nPJJ5rrFr3Hup7j/AIv/APA/5jpqv01n8Vvj935f9hWHEn9/L/8Ax/h/wf4f400001jq9NNNNKU0000pTTTXBh4KI/YPp76UrnTVr9SoXsZQAER+UBD5g5AewlEeOQ49xHj210hJKCiiuYvkAbssksAAdETewH6erpUDgeEx7G5HkxeO8OdmMgnPsD8fkc4IOKKHYZaN4v8A8g28e/c8fNXrXA+w/sHXiJe6wlbaprWSeioIDlVUIeXeMowy6TUoHdKpouXBBFJIpiCcyRlTEAxeS/MHMdD3xYNvtsA7bbFA5Z3lvm6yrOxjtmoxbqxx4/WMKMM2yE6m5iqjBkn103gRCrNGV9QjETCqhURaJEcWRWcbgrAfI5+cfyB81RmkGDHDLOp/FCu9c5xjvkHn1Ax6+mZOjmMBTAQhjEAQE5i9gMYQEBMmkP6xue4gBgH89YqyfmTF+DaDZcmZpvtbx9SqPDPJSyXG2vm0NHR0S0MiV/LqGWVMJY8h1G/KhBU7GL1dPbmPdnGeLRnsVFLNO7cNiDOCEiKTSoIPN5T/ACcyleoHLhSUnmGAj4tkKyi2BFsVm1uBZhacO4VPHDCJJyFxHwktpV6Ykabm43Ie9FGOSO2rAbub05zGnj5uv0Gl42gtXbGJJCxtiUQYqTqH6f1xoeIEwJ+kDrnYcEsQoHcsdo/j6jPI7/HBq64PlZlSTG4wucSheMtt5GBkZ5yPUV0MvGh8KR8gu+Yb89vUkY5TkScMrwCjQFUTFBc6TlFisRBdUTkFwQPMEBKQNXFTxaNstpEhtssfmHeywYj5drfbRseGy+1oS5u0PF3RQ0xX/gjy0ELIKwRgBz64kLKG6E/T/PtPg7attm2z0QuMMBYTxliShNJqVmkarUazHsIZvNvzIeslPIOmcUnzkEieaQpjgIlDg/YAGtyPnHBe3uDe2S+2Sr0uIbvWsa5kehoUU3Lgq5gF6SOKZ2YDggcVVDtRBEQKUR5U1hkmjjIG8NnkbTkYAyTntwP1/jjbstPv9RnW2sbG5url2Cw20UbNLcE42iIDIbcSAMHOOcVqCv4iWSrb1V3FXh+73z5Jl0zNaU2zTic+HcTqy5g/RDkTJiclcj0WvgACK8+WszwtzdJfQH6uQ88G4fxb0g5X8PbbMcSFUEqZN9r1ZUyR+A6hVHbgTk6fAAAdP/EAYTfoxIIGy+63iZFyk3XjtrG3+33l66MWRr+QMnqDjPCFxqYAbzJ2qZBZtbhKyCzsVGasOwcVFkD9sdyqo5ambgmraY7BW7XKzps4zrn1PHtYmzGm3uMsDtVaxaKZPrABo+vI5n9SLy4wTJNR6lILr0uBGcUIydKNmItwSPqi/W4JWAvEw/DKnh8DHPmJUjH9ks2GDFcZNeiHSlsFaXX7zSunHThbOW8MmoyugBZI7SGN5BIcNgTGCPerReKJCqHQ2y7mN3TLNtcwZb/E82I7Yc95GJCvIDaNJ7fVct2yjfikjleKqEbkccx0IcjnYlaLtW9iNTamM35RnAxMdx5OtqY7afuUykVRfdR4lF8O+iEFEaQOzaOS2ktW6EgJDS6l7ZK3HLY3GR62sf8AAFSrQvwNMZdPoffEuW3tl/Cc2T2fIlXyRccbuLteIKvT0PMzlrlxmlb+/nVohdza8inUaJK2e3s1IvmImjKNDMSSUuBETA+N05FHwtdiQHMr/dxxz5hl0XZzBEn5O5bFVKguP/Ed1EgWVAg88AChg+uinVl9bYAnKF2YttOCCcRsBgHGM/wOKztY9BKEUaz1NGyrtlE3T1nIGdfKXgdtbjYwsRuRiuSpAIBHOEB8OnDkwkvG5M3ibzc6Y7cgBbXijNe5prbcVZIil+7mv3uthVGJp+CXMVM7hkDxoJzFTHzC8Bz9l8IrwixOl5uzvbi6U6TlOqtFILLrGDjqcu1PiRfnKP6x+kQ5N799ZbJ4YuxQiqiwbcMfFUOoAHIWFP5pzq9QlWUJ6vpEDgUeB6vl4MHHzc660PDF2Lt1hIjtux0UnSm29Q0iFB6kRA4qJOii5L5CHJSiHSK3PA9g6Q5nfq39u0H3exkGTwD/AMI8fnzgepxUCy6DAx9fdQZ55/2cszwcY4+vRkqV8x43EjAGCa2RxtjrAuLanWqDj+vY4rVLqUGxr9cr0YnDljYKEi0wRjIxkkqVRQ6TdETEBQyvVwUvJe/OspMZGlRgKBHPqowBYxDLAyVi2gKin1AQVAbmTA4lAxukTciXqNx+sPMf8p4enh5wCBHthwhiCEZpNFVmq0qCDMFWyBPMerkB29QOduiQiZhUSIqcQOHUmQRKAx3rXjwUrQoSI29Ybg93mQETGevsS7Xcdu8hZLbV1uPlSFpWhJ6UpjFOrxLtZgylpQssdw0eSkWkkxcEdKKIN2q+rWhPHO6Tntn/AIXYc/Px2qDY9CFuNf18A85PTVkefn/fueceg9f1r9EA2au9w/EEHwIh2GRZ9+Pv+m+vv+WqNe010DKmVnIMpEvLOisMozOAmEBA5lE/NL5Yk7AQQMbq6x4APr+c0mynMeaZQIzHvhsbUNuWLrwoKlXzJlTIU1ac541rLzlxHTl62sfgCLr7m5t24JNpqhEzURmwerrJp2d0DUp1sp1b/Z7tvNoTa2fcBk/IljyOi8BNy4wT5+3HGLmIYLpqRbIMVsJm+tk3hiFMSblPxEY04JiqHatOjpGjyauB5FtHPH43Uemc5izxycAc57jFZItP6AY/a9R69Gvfjpe0diRjgBdewAeeSxx6rjipMNwO/rZhtXeVdjuF3F40xI9vBJmQqLC5WJrGvLE1g1mKE4vGkKK4PGzReTjfOIIpAHnpd/oOs8h4klsyuVKP2YbZbznN4qcLBVsq5XlUMCbWMhY/S5TLcMdZvSj8iqWNzMesjHtSihprIJ+IWkpD1rAY0EHOZZHwnvD/AJhU60vtjxbKKKA4Awv6+k74B2oms8KTzznEhXayKKrgpRAFDopGH9QNXBr4WexBiVNNntvxs0TR6fJSbQxkUkBKAgHkJkXAiQcGEOCABQDgADtq4bVRnzWfbH3pO/qf6nsfb2zyc8VNl0DwF1zqEH8RPTljg5xwANcBBHOSSc44ArX4lW8QvcOskvlzcLiPZ3jaeVJDWnCOE1y5Ly+1iWAGKpZKFuzPK0YlUsNjUMk6IQ2HJQK2kisx82U9b56F1Q8M3ZnOzUa/3GWW9b0ZKCcg7pE7u+yk1zDM45FJZJd4hQHQRleLAxs45Qj3ljQFN3615Dwp+sgNelTOP+662LcKB/d0x3+mIKaw/ClOVUzCAmIcRcDyUwgAiH1EAER10H8LHYgoCHXtxxyYWopC2MMSpygZEpipmTH1HymKU5gAe/Ye/fvpu1b+3Z+mOZPj/lfHf4FQLPoPjOu9QdiDjpyxGc+o/wB+nAweRznB8wyMbasH9BjXZloxeqNl1iqD6toeIZLLtzCBiIkWQADqolAA8wvJS9YJ/L7ceUyTnPDWIqJZskZLyNUqVSKPFGnbbZ7FMNGsNXoLzE0gVmHQKqC2QFVRIAOCanzF44D31r8r4XuxlbpFTbvjwRIoKpBCKOAlVHkOsvC/IG4EQAe/ACIfXVtlPCq2HScU6iXW3PHvonDVZt0pxQgZuVdMSGVbCdZQqa5B6VElDFOCapCHAoiXgY3at6NZk+26QZOO2fCOOc888HtQ2XQR/wDrvUQGMeXpux3jzDsfr0jITO3j7wB5HlrNeEty+3/c3R22ScD5fomWKNJTL2GZ2mlTjeSi1JyHP5Um2ZmDoOLxsooBHKYFEQE4dx78ZycpR7tIyK4M1kXCKp1yOiJnT8rqT4SWaqAdNwmcDfP1CAdRS/KOo+Me+E/ssxbTa9V8eUB9U5+rMGbaGyrDzJ2uWkJdk3BuW3Gt6TchlLe4TFYXc36EFFzuXJhRAVe1+kdomU6Mid5g3dbmhhaSHaNzq5qlFM4VNKI6FhcN21MeO6kihLuVStTt5wJU52SCLtEGa3rBMlUXOsQjNzHbgHssT7gAMYUs6xEn3J44BJ5IFY9I6YJMWmdU6y5YlVOv6Utkq5OAxa2vL8iMDYScFySyhOF3WTIuyKRqNwsWZNiVuq22TNN7k3c3ltF7Tj3PE+cQeLGkHDm80NvP1YHF/ZvAFnT8hEmOimx05aW/4cl/jgHZ12Kt9CkLda7gvenVKftU3K3t+CWNaMXJZr/TcxMyt13B1sX31xVKaWzy0AiCRL7FKQceFXeycO3B5I+s607URDxFcfuvLjAwnuUWkVTLqP5s7rb0NNVbclcJtGsXHZRC0hNHXKqo5WWhvh3oCFIk5B2YyHjck5hlrlQ7RhfepsTyFkJrZkTsrJA4orrPNmKpquPliPWSA22YWoEi4kFFGSK09EpwrZOOeItUyPnpTCctl1FZTgxStNwDFHG8jYOMEMqmNu/IVzjsTuDYluiNWkcJZahomtTyeaKHTdRtluHTAJ3Wt09tPEVIZd00cYfCsuUeNm3Rzlt/xruXqUXCWxSYhpmMkmlrqGS8azqELkvHdvjGrqOZ2nHtyTYOzxc2yjpaVikpX0zlMrCTeIA2H1AHT00LMb+No70G8pGjvt27xi5jEma+ZxXd1mMaTBAaPi4xar9FnS3SZGmCOmS09Y/i+JiirHyD4Iw3rgQb6G/3wDbRnDWG2uZfxzecTsXDKutNjGfEJbFeQcGNDmKM1GQmXYxlkJU7HHqrROoUjCH4RbRbBGSRZ/jsqcSRZ5Mbt03ubdt0TCfLi27t3Fspj2Ng8iY5sRmsXkDGFnk0XTg1OvsMk8eoxdkYHjXzWRZMX0qkg7ZLpFcqlKChtiOYSMECsrEE4bAIwASCM5BGexA7fnjk3/THU2lrv1HQNVsxwQZrVlBQkKZFYEq0YZgN6kqSy4JzWR8J7hcfbh6O1uWPHcqzXVePYGyVqxR6UVdaHboJYGtmp95rAvV1oSyVSRKrC2JgLlckdKgLUHC/Y+s7lETdPUQUymERT6u5Ux+gHTDjoPxyAh1CPv8AlrQ/NGw/AuZL1K5mghtWE9wku2Rh3+d8LTJ6NkiwQMOUBjanapRqRwFix8u+bxr+Xr6hWnxlWLjhVdoCn1awy1z5u22dSBz73BxxlPbiwb+bP7zsbxzmgzFWmJoSvGpclbemyFjYUzG1eIg6hX+SWGT7HIv5h3ApjT25JRZRhl3Rg4MsanIG1mwefjGT6dvcenNcUqyLukBiHvKCnPtgjOfyFS06axfjjLVDy5UYG+Y2uEJc6hYoiOscRPwLhJ7HvoWYbFeRSwD1JuWwvWyhFkSu2yC4kA3WiQSnAuQkVVVS9Y9RQMYeCnL0GAvPYOAEQHt7G7dXvwHtqUy+fKVx6nseAeD69+Pj27VUFWBYMuAcE57H2xjOfy4qu01RCuIHEgCbn7D7fbsPPIjz+X17a48xYoHEfcvPH1D2HjkP6+nH5iHDbQjHgHdwF5PoSRk9z29OO9VEkZWVw6sIfvgZJ/QYGfzqu1YJJTyknJ1DlRTIQvSoKgF+fp6uDcgAFKHAmEBH5gAe4e4dwO3InD6p9JRMcpeCFEpB8wBNz2Hq4EOCjwHI8D31obvAyJLWNOn7ZqBNuYnJOdZUrV3M11wdxO0LFMQPm3TJfw0QbA/ho9+pAU2SRF8y4NckFAX/AEYpKYblvBQeIHAdkTKhWwXOBu5wP1Pciujpdg+qXdvbRr4ltPG0t3PGGkFpaIoeaZgoyWiQFig8xCtxkVh+mYooG7vdRK7lch1GFtsFtqsstjPbg7no9NB1DXyGdKxmYrWrEqerTVOWdiWSFGm0XpSqwj6TEWwCuAklDKkC6x3HBQBRHy2xhP1gqQ5yKKHUT6S9AgJC9I9Q8gIjwGvH49p1Wx1SYeo06BQrFVqse1h4aFi2hGjRlGsSJtmSLFEpzCVsRFMhETGN1ikAdXI8jr2xzl4WVTIPInEgnbiB1R8o4FEAA3SBeeODBzxz25HnWK3jkgleaZ4miIVEKnJKgjazeUZPYA8kIB6A1sdRXf1hfQCEOLbRoYrTTIX2kNFbgLFIyq2PFYs01yBzJM7uxORVcgACUf8AGbuHcOA4+/cR9+/2/drkWwHMPmFDpH7GD+PHSPf6e/A/kGuj1ZgKUxkxBMw9JungTkMbuUVSDwCYdhKPAm6TiUvIh3DqKu6UMAIFIUphN1mcCCaqZCjxyVIvUCgG+gicnb5uPprcaSCYEPJEAe2JAcdvYDGM9hnjjkiuSiTswkkuSsgHIjAjjPYYKHcf1H8SRWJMy7f8S52gkIDJlOh7AxjpJCZjnLton6+JmWjdw1j5eMeB0qNpOPK6UOyXN5hUFRA3lGEA1qKrgPczgpd9JYIzRZczRSKYu5jGu5OzK22VtL1QxSNYyBySRs0LjmISaKu3Dxqap2j4mu3aAJmvT3kQO+RL5gqLAmJjmKmVTpS5Mn1BwIgY3JDiHUUwgAdIc8chwPlXt1qaHmkc2CBbLNCFK7aOJhgmgKyolFQROZcFR9MblI4nQIIKGKUSl5HjTkjtoQfDZGJwdjnKlsKQSR5sg5ZQRtznIILA+m0jXtYtkW0FvHqVizHOl6haS32nyeIQJMK6oUkmHkLx7HQ4eN1kSORIHJzK1a2xXSRyZSZ287ccoydnePLlsvy1MKw+27OslNP1XWTpbbyKTISxl3mLedGwJ5WGLkE7E5MvC/hSNC2A/iJZdum7ijbkYmwrViByBQr1TnUewu+JMp15Gp5SppZVFVzBSdkqhJKS8uNsTFueQgX5HpivY44OTJIiYqOvI5i3c7MavIyVCytk6hN5RNkklJwjvzJFQI2bYqKNwVUbMnKQN3zE5jmSTX6gAQA4lMABqEbcDiDBV2m67kLw3skbpKFkDHcW6ryE5tPgjWnHFdrcmZvYXmObXXJKwU8YllkKxwcDOyKqTuVF+WHVKKSPmiUNT6dLAzSTqDbgbVhicSSbvLhlUEELjcdpC9gOT5h6+36MTXVjnstG1XQzOVTwr7xLjQgzMButLyRIZLCI4RUtiLvxQxPjJsEbfqaZKioY/PQKvQkZ0AH6hTVMXkpAEA4MmIdQkN2ECgACXvzq5KcAJhEeAAO4/YOO4/uDX5rWXi17v8NYnopty208cWXyas1WoRbzk+cc40xbaJlxWZiUn7DJvo+u2takO3z2FUXi6ymznW5UXKqIzPW2TBzvVRc5eIHnCNjF6rTttdGp9sj46RjMu17MMplhzCRck0K/ayqWOnmPKalNi4IKLZJi8ssOfynRnZjlO39OefryxABzOGXsj206kFuFDFo1VST+JyqjuWArXv8A9kfUWjkz6hf9M2VtIW8OV+otLuFKqAzPGlncXMsoVWDCOGOWU/cRGfy1KSCiokUVVMIoFBUyZkB5EUzKlK35SH9YQSMBlDdQdgH5e4jrXrcHuOx7tyrMXP24bBYpuwzP4couN6HClseS8mWZRB1IFrdAqfr2R7DMBEx0pLmZ+sakGNjHzwy5RblSUiZzZlfxCEr3Kbe8KZ3j8y5eVTYxFnnqptajGGJcLK2hIrqODMWSksrvnlHfva6WTl6ySOq1jUfSMchHvCxormcIbHIeHRjemxNazDmvdfuNn8m4qxQ8q1pzg9ywvU2raLEYiTujuOiDtJRGqwk3LwbaVWjUpN+Zj5CDVN46AgHUuNSv5tv0a2C9jtJjOeQM5VnGTg+o/Cc4IrhnQ9Gs4QupdYWzhi4d9BtLy5a2wvka4FxBZNGjsVKFRIXQSY2lCDyljOV3XOn+at4VzumJsNRpwnartHQu6MBUH2O4cPjaEtuorqrRT8W2104j4W1qUsCxIYulmTiujM2v0/xNTJIb0Y/P0bEx2xMaJn9IqaTifvrezlgsbUCvv2J04yRCXbRUulM2gCOWUlHVHy4ostBJSi5pdh6fyzaFW8PCtllIqFxVhtPfHnG+VStXsKzjCGHI11slbuUQhJV/L+SJ90/hWqVJl1ZOKNYbamR7JIKWBs5PXlTqnQTznt88PnNsvS0KZleWqmynBykjNTKOzvZNPuG7FCwv5Nf44+t241WLrUpkmnZHaPZmVtGOnOLK+3YyMuk2RnXpIkq7vLF9ZHDXogLMD5I2byDgoSGQKxx3B3DAGPeoUdIaKomsbrV+o7tgGi+uLKHT7aDcfvOkF7dz3aAYMSeJaMuAXLAmOqk2ZcA7fcoTNTZ2bJW+Pe5+H2Kdlx5Q2yN+sWP4+zkZTASiFfKv5eKcGv7QeDUfKlkLMrWyO4gg/ERbj5uRUNvW8/dUp8U3bZYW23Y3cIlWS2zbXrmu7nhkGRiRL5O9biFY+KDI2O7pErSi8vjwMW1ozBWRZt/xE7+E+a938wlttwbtvpTLHOCsYVHF1GjnEo7Z1qpxibCOQczUgeVlVyl/SLio9kVDu1upYSmWMJgKAcAGY1GCKgj1l5AeB5H3AwcfMHbsPP5/UdbMdvBHJ9ICkz7cFjnBBxlRnPA9zzgYyQePN6lqN7qcjyXVxJLIU2xyJiFo1XaFhi2eWJewUoo2qCNpJyMDYFwJg/bfTkMaYDxdUsV0aFcSLppX6jEpR8am4mniknKqtj9Si6nqn6yi7op1ekrgwgAdgEM7goToKJC9jDx3AAEOREB5/P7/AJD9NfTeObNufKJwJhOYR+vKpxUU/wDmMPP5D7aqgTKAj27dhAPz+/5fu/frM0hyuEBC49cA9gfT0zkHggD1NcqNZWtikm2KfOd8XqOCN/A3SejP+Lk4ycD6AA4DsHsH0DX1p7aarWcDAA74AFNNNNKmrKexQScj8INLR3xX0wPRjAeNhkCshHj1gsvN9UDYB7Cv5Xl89urV2FZIOwnD/P8Ay/y1BvnPwcW2e/EuV8QUd12ecRgGCY/BylDwVMOsZWXyo525chKGyrHy7x36Jb1IlkKt+FfIeGK3OaUTFAANsC38MOLOn1f38fE6RMIiJyI7xbQiImDjqMuAQI+YuP8A94sPAn+XkA40pUo3npf84fz/ANNPPS/5w/n/AKai/wD92FFdh/v7+J/wPsP98y0d/wBn/gGuR8MGL4Ef7+/igcd+47y7Tx+8fgH8dKVJ76gvAj0qcgPHHSHIh/zAHV3KPfgefoOvgXrcPLHr5KrwBDAHJRE36gc8+5+/T9B4HkQ41E1OeHXTq02Baf8AEP8AEig0hSVVbjJb3Z1mkZFsQTrqpg7hm4rlRIICoVIpxJ1F6u5i8xs/iTbLcEwZ7Xt4vjPbzXUWqLKzNNqu4KyXFLFYLmMlDoZI+N/hUsMM2qg/LGOGQygvU4mVUVI39IkVwpX6iknrdUAADCVTgBFE4cKlEeeCmKAiAG7D2Awh+evNzt6qdbZ+rsE5GQiRyrCkEs/ZR51ztygddJArpwn56qJRKZQqQn6QMUBEBMGvzDxPhd+LBmpyCs14gG6LYZFQa6B49nWd3Nz3mPcks5Q4ncHmyz9UwEXGTyuIs0SALI1wCYVmlwOZgEMmL/bHPX+z9bSd1Y0odzed98+a3GOTviUt3kDc1NzaNedy3w004rXmrqDUCI+KrRMap1JKKCcWSRRD9AURdu/FRkcncAB3Ofu/J9sZHeto3Pisbcridw02yVnMW9J3HnO3uKe1rHieQW+PHKpjJw6WQ3U9O04YJCfWSkAiV2SUuKycRLHOml6dMrjxkKw8WbPCZ0J6d24bDG0SqDoBqLZ7vAkMtoSYj6laQbTrHAxcVu6+i2SACs17iM24mD+YaPCDSGQvLfwuYGMQKEXvo8SRi1bogidshu6sbPhFIpU2hTrpQA9BGxCmKAdB+op/8PT3xFlTahhHCUGaxZO8THxDKnEElWsWou43oWY5yyb0FTJx65GkC4dnduAQUEiijZMC9BgP0iYoDWSVIRyN+0dwfujAPOcgDOec+hrZt9MvNRkjXSZpNTaQhQgjkfc+QCihVYkjgcKAQR2PFZxQ8JralZ1jSO4lLI292VYKAtWF94l5XzYjjqRKYqsv/ZylNx6H4NJPLEYnmUWhnIPSQ8QQ5xCPIY0j5YSNjSJIMGbZkkm3TaFRQTSRQKggQE0RFuQhSqFRKAFTAxyiUvIAOvzvSGDb3k1RVntbzZ4sV8kXJxka3eMi7xrVjrBlqrSfAKykDf2dXtkmqq7EyBodqrUUQfplcnUWamQAqubcGeGJurjbg3veavEN3lvKzLQsk0fYLTzvY7K2q8q/FqZp8Oym7SjXVkLCgkumi8WpcQeSByY5kWgo9Kmo9/dSFRblFU43AhlcdgTyFVgccBGY4AIXA49UnR89pH42t6zpnTcnGYY7hLjUpFIUhfoNuJHXdk7TcCCMspQzLJhBK9lTP+FMLQg2fJd+q9Uj2cknAjIPXhQcJSsiVRQI86DIrt02XXBicxiikJS+QAHEo9POsEnvIuWQDFhtq23+2X+afOvisVe8nnNjPDFnrLXqI5n67kBi0uUnKPlzOmakJHr1RkEm1UeLqumRmpU1tcoXwQtt0HkuxZhr2eN71Zyhc3My/stzq+5ibYScg4n127mVUePS1/zTDJuG6Kq4CU3mGQKJh+QOcut/C8hmZEhDfB4mKYplIJSN94VnU8gwc8lIAV9MPL7jyHygbgOxeONUjguZ1P04kRgll2Oo3ENlC2Q3JHJXbw2DuYZFTNe9C6SQ0FpqXUN35VF/fIdNs1l2YbxNPt2ne4h8RQfFN8m+IkeFG2GHt2237dPmUSvM67gHNBq0uX42OO8FN3FVtFOmTAItYVtmtGQLIWqCjU1XSDsFqVC/Gji0dKJMjNipKZbxLsu2/wCH7SlkGu0OGlcnmh30LK5dsqCU1k2woyR2q0o5stnWRI6mn0yszbOZV2sVMzhdumoYvIdtc0/DIbck9Tv08TU4gmQUVP731oScl8wOViGEIQ/SAiVP9GAjxx7m+tcXwwY7q6v7+fieAUQ9g3lWkDcfmHwDnjvwId/562ltbcIBI5K91XcQMjHAyGxnOfIRjPGBiufedXdQXcclmkkFhpDH/wAlZwpZ20qgAIY/B2s+FA3GUO74DSFm3EyTJxbNJIpG5EUhAhi9IFKKZPbp5TApQMIcCAdygXn27a7Ei8nEBU5X8somTBTzgKAgIcAUSkA5vupyBjdx6Q41Gz/uw4sf1d+fie8CPfp3k2kAEfz/APh8OP8Ap+4NdZ/DHikzl/8At2+JubpDq/SbxbOoPbt1JAECXhUvPAKc/JyPYerU7IGPlXnjBPc8DAPbj2BOfjnNebklFzzcSsSOc7TLgkg5D8BT2yccZqR5nIwbWVWjAk434sRsRRWMK9bA+K3Q5AFysPOBwRIOv9IsJPKDkvzcca1szzv72b7ZF681zxuEx3jdW2oyStcCYlFHhZYkSZoSS8g0K3lPLFmZ+zKqDgEREVyeX18H6Yk2HgHw5d+cnvCd78N6z9k9w4ji1WDl8jyq+Vnibldsu6SWzyMmaRe1Nymx8h/RzU1NscqiPMsby+kZb9vOynaltMJaG23Lb1i3DJrstFGta9BqcfCmsy8QR8WMczZUSB5qrMJB8LY5jHEovHPHHWPN1UqMH9PgHsP5/wBBWP7v3y/OMs2849BkgYwPTHFaqye/7LWYkDR+xraTkrLqzjpmatmHNbg2B9q+Qcek5D8SUTM0WwyVPTDuZ85i5qDJ3jmN+MMPibhw4jTtCIOPuLxD4l+4tFGSzbmug7PaLOD+HrhgXb23d5KyE1hGfHXaKPu1eK43m6faJwT8pooYmfBWytjEReSYPDGRk6YJqpmXZA2AiCiplC+YUE25kVRHqUbkIAgbpEpeprwQqACBSqKdYiF9YkUTSOVZTzDgqoIqcjwcO3zAUf8AywH38vkwE9uodWq9RktPCR2oWJ0lIbiWmRt6cjErFc0qQ3mX1xnt7jVUyqa8iXHTidjWJqwjYV28c4sCbYFAkl4WHOoICxT5kiLBNUQ8tuwjkkhJ5IgRogQfSkAoEa/IUo9JhKUxjdXAeUUOkeeQ9DyH3/oPf+GmlRke4/eKsyTV6l1/KiYpAMVMnV0iqPPCapT9I+lBIgn6UgBUOTh8wdPe4NW4tkQRFU63SY3Cio9SpiiPYVTj3UU/51B4E49+A9tVOmlTTTTTSlNNNNKU18n4EogPtxr618m/VH+vroO4x39Kfrj59vn9KtxU0yiIlJwPPuI8/nyHIfmIdh/66+FEFzqCciolIcnSJB+YoG5DgxQEfl7cgYOB6uQ7hwGqgoD7AYfuPPYR+/t2H+X8td/IAUeff6D9gAR9vr786tKygKZBu7Dnnttz6H+Hp29qwtF453NdNcKDnah24weO37vj1qwqMlCuAFuYoJAkomQqpPOMTzTFNykInL5aafRwKQdQKdRRE5RJwPeLFYwJ+YYFDlA5DiU4kIoURASqClwIdY8cAUTcEAxg6h1XqkBUOP1uDAIBz08e/fuA9/p2+muFusBApRKADwA8hybgB+g8gH0+/wBvvqAts+1UVFP4vLg+mMnvk544457jNXJJTwmVbaBSCrmZhISMZLALuXBwUAY7ue2AK8i/qUJJtnzaUh4x0hIorN3SKzZFUh0XKKiLghimT+YipFDEVDn5wN9B1oy78L3Z8xYybrGeJKlhbITlqurC5TxfBtK3d6zPnbqpt7HCSzQE3DSRQUVP5hiqda7RZ0080gLmULIqREEw4J3AR5HqEeePyHjnkf5e3P11wcRMYSj2DgR7BwP24Hn+f9AOGe1jbkkMoIx4kPj7ckcqhPlYY4dWG3HHuetp2v67pAP1TrGqAyACQx3bokygcLOnC7QCQocMpBIK8mvzVucb+IHsIcrzOYd1VnveG0PKYTm8W0vJW8QVQJKgL4Au+05dy2iaDUIQ7L4MrluOy7ZZMr51GsS0vy55deN3ppmSt8UzUYa64qtG3DehRb7DRtmirWlNOMCIVyIfN03MQk3YxsNlFezpT7RyD31kmWuuYz0RUQYuRdqHbyvPmhFUAbm+VNc5SH5OHIh3Nx0iQxVAESAAkHgBDkRHkADUb2Q9iT6j3ezZq2M2Wp7W8x3uelbFlcEqOjasd51dS7haYkX2RqYnMVtBe+qSiZGVeyIq+fPKrFTNlRQh5EsqqUmr9Vh/NFLJlQCGE7wYAIyNjCZGyOMOh4PJBw1d5OuL2T/1fTumtQLMGlOp2UbTsRtAK39ssF5EwOGWOKSON3GJFZHdGj+ueNrfjSw2fIOPcDbhtpW5+0WKVvVzndncQjmXbPla1zDpZwwuGfKu9k8RhmycaNH8sBomULBEg3ss5FGReAIifOu2DxeMbDZEtt+6BxO1DNtQYum0zkQaStCYyt4sHzKLjRf9L5yNFyhNlci+suLGw2OLpj1F/ENbtYEWJX7jYzG2+D4NaqvhLe5V6ZtZ3O2uUBKjVGMvj+/Y1yXEu0F3vqsZZLfVKnHsTqBbpIN7swlK1X0oGUfMGTZ5LFVFdPfU8ZXbPHrJKtIqajZePc9JxM3cNJGPdFKQyxFk01CLouG6phMoUxgEFAKUTgcTBMi6qiqUVVGVLcq25cYICgKo9CTk/pniV1bo+7ZjqHRF8H8M+Hqmk6+bBYGLDBNtPY37Sxpgn6MZI2bygy+UtVQzvVKeOE/RWSFcrueCJIIyzFZZUTcmAE0CLmUOIAHcClEQ+3HOrkS2wDpZ01j5aPeuGKvpnyLR22dKMnBgN0oO0kVTmbrD5anCSwEMPQcOOSjqD3cL4PON63bnu5nZdBBUs/xTKNZwuPFrg6oeOZmMRR8l4jTJ6Nh5pfD9+fmK2eq5PjoC1PSN28pE/ATpzqjpnoZsn214ywTnPOtO3Qbid1e2bcDuly2jP16kmnHmMYK/WCNLY17W7pGXWcrLq7i67CvZJNuXKdmq+LX79N80eDVG55tVuyotxfLxNDEzLnJfUPAOw4CHZ9HYjnABBPuTyRWaSy/Z8Y47qPqrUmVwI2ju9E2zRStIi+C3g306+FGhLK4O5z5SiDmv0/2XItPq8TPyj+0whSQkc/kXqBJJgL1IGDdRw6KDE7kgisQEjkKgc5BMcSkHgfaIrA+9Pa0lYpLcJnvK9Za5QuYPXeMDSir1welbf7o6Sm8eQzaPQZOUalPXCDaQ0peoxm9fpPpuGbKLOBM1TMfw182L7YMpZ6ZYErOK6lZW568W+7htwB2DW0ZDVWWdRq1Srdlt4+mWSumW2r6StadxAztRdjWZtsaOWLKGWbzEYvw5hvANDRpGOqfVseUWHXkH5odhHMIeBRcv1QUdSCgmMRBD1K4lMZwYQ6jHAekOe/MS31K6v3huRbtb7PFUxanNctlWXarRC2gGFyQXMhO5sANtJr03gdC9G6E1vBJr2u6jrUySyW1sLLQvomksrCWV75m1WR45nkRjZm0jMsEau0qiVFbWxHxCMPTIvXeJablzcIgzXFCxSGDKalbW9cfkExGsXZ1ZSWrh2j10QrpVi3RI6IdJk6EypBTAD2GW3V52yYdGF2+bUsloWhsB5KSc7jWieGaS4rhRBJy2jbVDlyI6dWMJBeOOhFHhUUlmab90L5MzUiK+zmQc24QxBQLHkS/ZJoFUp1VYt5WyTa03HKNIxoKqbczlVzHHcvlEyLOCI+YhHKq9ShSgkBTnMTT4PF42AzC6UPhzNSOe8kSAiFVwzhWLkrHku+ySQCY8JVYWRRgmMlKJsQdyijNzKM0hbMHC/qBOkRNTsG1k2bFyGUgLkq0WPUhH3eYD0LH1JB7V5JdY6Ss51udP6chlulybNr/Ubi5t3BxuWeKFLRZMnkPGYGyODnzVe29o8S2yqmr0hiXAGKmsoRFEMhsMpy2UH9dU6AWVckob3H1PbWAiBkxj0kVrLF+oScetE6RkfTnvSeG9/kiiqwlN2+PGEa85K6ND7bWUVJJNVh6XKLN6TJqpmb8zY6oNX5SKnZufLclRVFICD4h/4g+Ybe0ewOGfDx3lssmzSZW9UW3EY7iMLYbGRTORQimQMqQ9nyLJ1CLNHEdlQfNKTPnUkTMmItSEcnco+ZbX3xhchuE6k923bTdurWbP6dTNcbuMsef5XHRQKK4T7bD0nhTGTG8PAKkMUWEcXuupJIyC0gEiczErRzYaZblcttZsgrJmRCv3c4ETxxnnJy6Mcj73aqydYyOQ9r0f0zZBQ4a8i0a2nud5xl1m1A3lzGQAuFhmjjwFIjDMxbNhtg0Qoh1n3Qby1TlRREVFM/TSgKFKmHYTGjQFYwgAmVUMBB8wAP0iPtSn8NnZA1RLZbvh6hWabBIi89kS9sWMhY7E4XEp3ExbrE9BI0tJyLoSvJSSdAgZ9InFwcpTnANa/wAhhjeNIBKweePFOx4XFaSD9zkBDEmAIHAeWixMYkq8cHp+VY7NdtlaTLNVm6Z3cmhWpJyLIjyN9H/xplE9FU7x4WN5H8JJb9N5HiAKTZytR2cr56t2bkc9JFKLkao4xTZK/WIK9gyMiWfJEyc/FINxhgflcmVZJoq2XT40AMlxI4zxuKtxkDgOWwBn0HA/QHG37QOr5FMdnr0y2ykDdb3H0W6Xja0RihijLRrgBA0jBgd2AeKmzbyezHaJiiclW81iHFOJqUk6sE24j3EQlCwwzckgL6SVZRhnjs60nKvEFFAQbLHUWW8w3SACOtcx8XDw+ZZwhEYgy5HZ7yI8MKNYw3hCIfWLKV7kS/OvHU6sSLWAj5WebtiOHgtnc2wTLGNpBYrkTpFSV0DpOFsGQN5q952H+BYTGWd6DJLWdped2WO4PaxVoSLMzdRbt1T8i1ZDM0i/vS60m2RbVw9aYN3kQtNPDzSB2CbZ5uvC1TxV82LfGpWd2xbC14QnoU4CrVQ28xPIplwBQZuSnJ1Pb+tj99AHQKwJDNI+zpyZX6rk8iz9GVFxY29g2U5kcBcbskBcDjGcYGR2PYeuMVwZNUv7uR3uZZ55JNzNK0zy3BkY5JePIwjZJaQv944KksxGOtwu9uMzbiq94Xc+HjujLLXeFcegabqKHWMGYTTWiV0Z0svecns7VkaWprWDLGjLRs9GUSzOELEziW5GIC49Qh+bjHVL8VgIqcs2HsLbjdsG3e0bg4GZte5PF2WZvcFlPIFAs1UvVht2cRQla1iNTOWN5OeTipuEvczK0qRx9BSjWNbV6dGSUVbz1ZS8EW2bod121jd3u73aG3BT23mQsH4uxM4xChXcN5IZLLkGvRkFSD3+bZ45PBuGzCVnnLdSy/imdjGssqlHKJlS1P8AlZt3bdVkYiZm6XQiDQUClQKgiHSmCSHUJE0RApegS8gcvBhITjpGVsrQHCqF5wCD5gGHbgggE9+c4zk81s6br2taNKWOoLbROfPFE30jxQNoPiQSqY5AEwWVsgHGCGXK/mno2K95mJ6xXoyq7lk7z4emRolhbJPK203FqTTdZkKzXJqna7flSw5DC6MndHkLLMi8slyycBLBNPnLtw1cVoh5tdVnmRHZbsg3r42yBi994hO+PMuMrA1YV/JVZs28SxuK0/buVUpdjHSzSarjNN/+njE3DYgkMQAbicwkOAFHdnIO0634JnrVn7Z9P2JrJyc49td92qS9nct9v2TyzTxSVv69YrCDNdpQMt2mZU/F62Q28fPuZuVaPoBxFtU7SvLRmnu5+h4j8WygKYNjNvlox5kl9HTMRbbPn2iw9dv+3OrmeoLJz7OsR85NllrTYnyEWpXoBOdi276vLyss7mY9zHpRj7Ula7tZtsUQIJZg+9yi9ii+beUHoM5wOCUUFq9dp76F1Y8jajpiaV4aIZtf0oR2+krEgxFdXmjhFjcxtmOUW1xCfElUQwzu0dvUrW0Ta5gLaNhSk4M24VKPp2N6lDpNYsGItnD+bWOVEzyzzMy2ST+NTdkXJ8VmZcSEGTeuFnRiFMrwGzXwxRRwissKJvIFRZHlMBKk6MfpK4TL1fIoLcyqZhAfdUw8/QdTtiu2ENm227FW3BLIt7yujjOpxNfC7ZBlV5KUkRj2TZmKcQ2XVcjX6y3FHy4KrpvX6FejwQi0XzsiALH3M1viSWRVeYbZCoyP0r50yW0c10lqsqxLcygmVWXe4I3SRqxO2N+CoBI74OKaaaamppppppSmmmmlKaaaaUq1GEQA3SQTFKIGOJewKG44MZNIexje3PJg+n568hNW+ArDRFeyWGJiDonVTTXmpFnEDIgkTrXWZkXXHzm4l4HoTMbrAA6hJ8vVGfnyjeLbe94YscFZq28YX2PBjCHUVsc5jActZqHKiaz0s1GM4F1O05mjAu2ijYWk3+IDLx6yQmLGL+Z8vrYrwjNpNjSZqbj0cjb1nEOZJWn/AN8C7LZtZ40cnMCsufGDObYoFpqM+qlHnmUGii4PCw8ORQ/DAgmUq3/71rbpbFVYzbFWs071XcKsLe2/3VscJ3iPx0usY6UGleFZ2fqXwRGdWRkSxrpknKiqnESih0k/TplXsDZDxZc8pHQnZTbrsSRiBBYjSqC93gSGTI6WE3W4efHozAxcbP6wi1KAtWhrWWdcTRhUcR3wNMZCV9CGjWDQjOOaIxyCaCbVIrBJJqKaKROhIhBSIAgCZAApOeekvIfUdVR2iJ0+k3WKfSPUn1fIcPqUxePmKI9hDnvpUNgrgnj17Yxx/n25+Ki5T8KPalZ3DZzuYRyhvVeMDC6qbreLdFs4ssZu+pNSZPjdObZthqKVlWTjzyiDX1BZAkHFEUOAR6ZjSNRbKPjW4oRzQGLVHy26RUmqLdqoiUoESMfoEvnFSIQAKsYCiQpuAIHPGsP7nSZ8U295aDa2NUZbghp0uXD6l7RO6qLO4AiJoZexM00znUYAqU5XDcgG6+pMoH4ERDQPb3i7fnn3C+JJXd5n+wUGRtVSi7Dfce4qxwXAWSqldlimBGCeWuLuNrV+FIiDn4nEfDypvw9Ic6qPkgBqSifw/wCilVkyNzORtCkgHGSMtyMAAnk+1d3S9MtLsFZtQtrBER5JJZ97OyIAxSCONHLzsAfDV3hRiNniAsAZAMqbiMN4XhXU9lC/V6tQzaVCHWdun4OBSeveoWzNVBqRZ008wG6vKPkmIQS8GUKIhzrPIbu8jZQTOy2o4NueQHMk3CYrmQMkNP7OsGWevs+QO4isktC2uWMvI+ekaBSUpqfq00np1DNxSKCmasUbLtveJbMnf4SgxsplN1GO4+dytYSozeQrP64yJpN1aLSsim6l38qdFAzxwsiUy50SibuAc7Qtolm3D9C3IgUAMAkSApUjlH6GIBQKPSIfIH+Hkfv21zFdzKGLruHIZiMEn1ZCCD6gDbwfNk10vrPo3SJPDsdOudYvAylptVzFpU23bujn0y3JlKsy5E31gpZDgRRsN1R5MMIbvMtk6s8Z4Vx9UZ1MJlTGmAmy9XtdMlSf+RCtM4JP0n9ugWpVVCPBcUmC+Kj6dQ6DXyQIOVcU7KNvWJLYhfa3RI+Xyidg6Zz2WrKuSXyRcSSxkjTEpbbGs1K4mpaSO0bGdu1uk7g6YCbp6e+4IJoKFFPgekwCUwHAvSYo8c9Xvz7fX29/fgNdCqYNwMAqgKYAAmMoHmdHAcGUUHkOevsAmHjjp+v0pHYCOQTTMI1UgkDcAW4JbLM23nnCkBeQoC8DUm6v1W/hkgsr200G2lVo5LDQo4rG0liK7fBlSNA0qlfKxmMskuxTLJJIN5+YsjcixEWqZ0CN0CCZE5QL8igm8kC8cAHT0Kcl4EO4dwHsPoeA+wfwDVkYCYzk/JRIQpEykIVMATAS8icyawGEVSCHTxyQgl/MDfLfNbu+CTBgZWRQFypBwQORkevOf1xk4zXmV3+YyhN7MSxVmYOSBliW/ExySB5fbiuOA+wfwDTpL9g/gGudNTV6+RIUfcoa56S/YP4a500oee/OO2a44D7B/ANOAH3AB/drnTSmB7VxwH2D+AacB9g/hrnXA88Dx78Dx+3SmAOwxTgPsH8A1wJAH8v2cB/01QesUAFBOkIGTABFNIQUOPPPBAA3R+kHgeCgI+w/UNVKK5VCc8iIgIlNyUCiUwAHUUxeexgH3ABEAH66A45FVZdwx8g129Afcf5f6a5AgAPPI9v6+2qN68K0bmXHsUgh1GMHyEDgfmUEORKQO3JgAwhyHYdUqT1cpQ88Sic/kETECgVBVRQDmEE1e4iIgXkxRIHRwAAJurkJ3H3P8/z/AJ+5rGIe3A7j0Pfj+P8APqavGmmmorNTTTTSlNNNNKU0000pXHAfYP4BpwA+4AP7tc6aHnvz+dQAB2AH5cf5VxwH2D+AaCAD7gA8a500HHbj8uKkgHuM/nz27furjpD7B/DXHQX7fzH/AF19aaUHHbj8uK44D24Dt9+/+evhQpOg4iBQ4KYeRDkA4AR5EPqAe4h9Q7a7NcG7FMP5D/lpUEA9wD+YFYZyhiXHObce2fGeV6XX8jY6t8OWKsVUs8K3m4Oyxyzhu9K1cR7kxSKtU1miKhWpjAQpyJm6xFMOY/ZHCO8TaIVubZy4oma9vEUKbCu7R8nu3OOprFtfclF7JKUHNbFndnT6t1dFilBY9wmOP4iLYMn7VmS4JJxhDOZVDCboFQA8wOkxyFWU8sAE5imInwBT8cF6uDf4fYAEBEQ+0EEHKJlDpgBlhL54AbqExkw4BM5uA6gTERAoccB20yff+R2q2TjGTj29P3Vqdt93jYL3Ht5qPx/bCM7zSHkdAZJxJaCEiMk40trtB2oajXSBXWMRhbI1Vg/bS7CKezDdq7YOU/WqFIQ6lv3kSOGUNumSCZ9x3C5QxzJxJIiy0KXim002tBH8kyQi4VWNcgAP2qkoLJdwYo+Ykigd4CCgoeUN03H7MsGbhn9cvFrqzaLzVjxjJJYfzvAopMsp4hkpNZq6dy1GsyZSvIpV+4Ytk5pNFQgykeLhioqkVwKhYIMd4z8YbdnnC/03LW6bbVX8VbXc7KzmIc30Pb+nPJXa3VdO01kIueoC12h2jF/GMJOQLdYglmkE4CdBmwRdShDi7JrzsrgRFXYKdzGMAsowSCckYXvnzLu+53INeg6X0u11HU2udSkSDQtFt5L/AFcuTH9J2tGbWC3kCsRcSSYiQ4baGeTY/hlTQYC8K3xH6Tvhu1vJudU2+7FMlYbq0FK4pwzf5obZUpGHbxYQWJsWOCx8M3psVjdiWRgadlRum5fJV1s6iyUxmnYVzMJIpXYFtsTbyUfmveBuszhi5sydOb9iXcJuxc33FNyimBBWcV+947koxJCfr5DplfrxTpVMSP2DA4EMdIA16GM2C5jyIgNsz/v83WSeRHPmtVpDbDkaX2w4rPApKFCCVjsSRLy7M4ycWYlA05LEsC6ks8/4hRFAfk17GL8Hbw3DPmFnt21LF2UcjIvm83O5byhXmFvylerYkqR04vV8tzxuk7slzlZADyktPOE01n8iss6OQDKCGrWyFIs8ZZi2Rx3Oce/fnAOBjjsMc/VdTudYv5tTuyDPMfAzGMIsVuqoqJgD7LZsEffIHOSM1HJjq0f7PVjm7Vy/7WcQber7n+qyYPsbVDbniBKWzlJ2BRBdBKOx8yfEgWik2nHKvuU1JdoQWZHQ+aAgBT7svt5+5O/LJxO1bw5MzRV6YpfEH5t6DVlthxu3q7fpaLfhi71hpmd9IW5R65jioVdSuMGzmPGUkDTSR4xJq7llQg4lsdI7aPaNzIAIJ+S2QT6QHtwAkTAQAAAADgQ41Wg1TKPJBOT5hOPSbjqEfobt3KHPYPYNZq5uBxwOO3Hb8vbsKiOYRHi253K5eT1z22bFDwaqbRKv1yAc71WuRAcEOu4fLy1hLt+PR3FdWQIwGMbsLEnOndi+UeRxmZW69/aeH7lS+AjY8+b/AHd9N5PegDeYebcr/K7YcXtIwgf+EJxOJYiTvTGAkiMyJpzMijYXh5p6Kj1RFqY/lBKmo1QVAQOQDdRBIIj79BhKYxefsYSlE3Pvxr4BkgUoFKBidJQImJTdJkiBxwRMQD5S8FAOPbgONKtk+5qLJp4Q/huEfMrNYdqOI8g35rLtpady5kKJZ2TJ89aY94m6eXyyW1wwTdzl3np1MkrY5hZNFaRl3Tl6ces4gMmraMYIKEWJHtUFCdXWqKCPCJxN3IgYqRTATnsUwCUQL2AAAeNXsWqJhMJi89fT188fOJeOkTdu4l47D9OR19AgQBAQE3UA9zc/Mb8jDx3D68ffjSowPau3pKPuUB/cGrc5MKQicnfvwBQDtx9+OPcPp/273LVtVDlQQKJQUKADwcOwB9eB7/Nzxx7cd++skQXeCQMDnHYHkDk44HufiqO0SbXlDPsOUiXlpXGMIFyM57/ABNdBDHP0qdIdRhEekQ4Ape4GH29xEeP+vbtUB5JSGLyBS/MYerjgQ+vP5cD3Djt3+w66F1vISBUBN8huk4GKBusBKIce/BA6uB6g54DnsPPGtOMp5dtl7tbzAuEFvLuTdFsvkXJSDYH1exVX5BMDkb9RFkfiGQZhFVMa9A+cyTTjyS06vJouoVKMfRLOPMCNq5K5yDgdwcEjO3j8s57A1v6bpk+qXH0sOqQBC7yPnZYxqVVm3Y5jBYLLxulcoFQnAPVl3K1yu1yVwJgI5Ura3TaL5Jyaq3B7V8SVd0kVUEBRIoQJi/zqJ0yV2tC4jUSxwythVl0XMEjFSMZ3iH7cYY+OoTCGxi3X7GviN2Ot2ZTCmRcZ2dzVLhDwM1YoicyVc9wt5jElpImJrJOIxi9senZu3M1kWQpb34YYvmLtt4bjeYbbLAwO2zbnVS5E3H5FJLTtXrruXM5FMij9FK2Z8zjaSR6zljBpyEkk9sE4eLfvrDdpqDiQaItppxMR2we3PbnC4NhbBOTtjc5KzVkhzHy2aMwy7MGsze7GwQWSbtWjMzp+pAUmtg5dxtEphZKSb1GvC2hWz94k2KufWgtyGVFDMN3LeYggldpBIz7MMjIIwTjAHT1LVoobf6p0ryxhVBT11BipBeUA7WBywhiVtltGWwWeR3ewbCaHuRxfthxPQt3GU2mZdwFcqcMyyDemkeq0B9MJsGxHDRZ+u9euLQ4YLkUbq3FyWOdWkxTTC8PFKOjM090NWtmmHUBiAIEKHSUBHqHgRARER4D34D9/Pt9bprak4bHsB6/6nHyO3tXn1WVQBMwaQgEqDu8PIGIt/wCPZ23EAn1ApppprHU00000pTTTTSlNNNNKVbvRiJTAJUxOBwEipw8zgxeeFxIPTyp3HkvV/wDqD21UtkTIEOUxynMdQyhjESBIBE3HPJQMbk3bubnke3btqo00pXWcBEe3Ih/17/T/AK6CURIBfqPPP7B5/wC3bXZpoOCD7VDAMCD2PtVoKxXKY5yCCZjGIkJerrRFuTq4OCQgUCKm6h6g6jAHbgR9wAxOkfqR/RkBMpATA3UmHAjyAJiHHt/i57+wh2Dm76aHn1Izxxx/PzVGiV0CMWIVldSGIZXXsykdiOcfBIqxkjhQUA6SZREU00j9xKQASEwlEhAAQKJus3Vwb2Av151cQSMA9Yh3/wAQ8/u9hD29vtqr18n56R49+NCW243MAo9D/Pt27ZrJgZJUKjuRvkVQHk57yN+L+B+apDdAgIfMbn/m44+v1AR/y49tQ0bzN7m/DBuc7Jj7B2wGZzhiuKodZssZmNnc1IdpL2KTNMBYaqEcSAfenUrqTGOUUdFdKC4CSKHkkFIOqZYpBEpimMIhzwH3AQ9/r+erDINCL9KS5zHREwAJBDrKIj7EUTH5VET8/pEzCAD0l5Ae3GndQ3OpWUkdpL9HkOdpdVbO0c/fjdcn1DIefUCux03qtjoev2Wo6x0tZ9ZWELhpbC5e5ghlHA2u9ncW1yu3BIMc6Ht8Y/MDsX8afcrlTcXb6Tlvbza56lPI19MwFVxNVwf3PHBY5YCKspsJp3XW9ijn/qUyDMrOIxZAzA3po9356ot5okPEgwzDoEXy9T8ubdwdgYYdLNtLRglZ8UOPXliC1eauALGjfMa+u9SZoCfrG3lCt1KeXtNFY8pMRdpq8QtKgmVrtEdDQ9ktCTRFvOykdWRfjBR7lVNATu2kUaWkjR6SzkCIGeOBSAvmG590qz5IYCNygCZSkKoqoKS5u4+YYyhSKCQC9uBATeZ1d+npDniaFpXUOl2bW1/d2t3P9ImlWQRRxJ4UjAohSIhiVBOPMMfdHkUCvrv7Tuvv2U9Y9RQ6r0z+yEdAacNH02ym0XS+qLqVDe2sAjnv45LuxlRWuZB9rmB/EKtK5E0sjjST/efbLf8A8aoL/wDarR//AEWtisZ7msHZcrRbbSMkVeTgjvXMcDtaQJFH9YzEgOUvSzAMHgAn5hP0gtwTPz8hzcDxmXyEe36Bv9P8IfX7/L9PrrVjLWyHannW3KXvLeCMaX23rMGcYpYLNUYiZlFI9gCgMWRnrxA64t2wKqAgkJuhIFDAUA5HXZYaooyEsW7cO8sI9Pxgz8j22c88jjPzpZeh7wiJ7PXtGUeZryK+s9dk4AAiFibPSBhyeZjdnw9p+yffldgj5Nx0mUTGvVR4D7WGJEf4er76vbSzwT9si8ZSTV21cAAouG66KyCoD3AU1UlDJqgYO5TJmOUwdwEdaIqeF3sI6eQ2rYVT479RcfV4pg7CHBRBr8pvsb3+3YQ1aDeHrjRkdJKrZT3LY6rDFuklV8f4+zfZarSqo0JyUkVV6TGlJEQDCKJ0IRkeyW8puidRMgkKXg1BJqmeYNObgYCXz5JOOMvbqMn2GT+gJq5s+gGAWPqLqNH3YY3HTcATbgY2i31WYltx5ztAU5GTwZEiyrI/HQsU/cQ7GLzwHubjnuXj6+/5a7ivUDlA5eoxBIJwMAFEOkOOB7G5Dq55LyHcAHnjjUeSmw6vJFVPE7j94rSQIQwILr7gbe6QAQ+YpgRMqQDplUKmJkRMUpk/MREQKoYQ83JYo301tNw+U3q0RpGMEk1HshJbcKsgkdNv1FMo9fOsiJJM24gcOgpznTKPUJjBz3xfTtRj/r9HkODyIbqJ2wAMkeIYlOecc/8AQGE6c0G6Yiz6x060crlY+oLS/wBMLkkDbALO21MyEDDMZDEADxnzbZMjPESlAwG5KJujn6dX29w7j9B9vfvqifTkZFsnEhKPW0axakBR29fros2jNIwCILO3LlRNBBIOkepRRQCh25HvqAhXdzufLlN1j7DOe227ay021OK9kuMqG3qt0/DlAlEjnTYQ2ac0OL+6f43hJgW8iZtZ6RTclG/8KdcseASFXRTcHkHe3kbKFipG557hHcHASjl3K0LZNiPPl9iafdqjYFRVcY1ylAMsTO6huMjXANI8sUhdXFYTYA2fGM1IMkoVPXbX4lYK9nc2xOcrcKDIAPxEQGZSGI2oAd5JGUAwx6+l/sz1LWZWttL13pe+mwpWRdTksoCWK+Qvqttp+3Yrb5Jm/oyICfGLgxiVLeB4oNrpTWZqHh/7ZLb4heZ6PkFpAZbquNliw1RxK0Ko8SfSdptztJUriwyRm65qpGQUbMRssnHTJ3c3GelbA92qwl4kO1HM9kQxkwyAvTM2Mq+6m7lhTJMW5q+Q8duIszRGarl0IczyttrFDu3iTR4zi7HLJnV6hbuFyEMcI4cW5d39QO3+p1+nbOYTZ5QqOZtSqTTMW12LzDkCqw0GiKEfEROGplDGdTaUt2hwCc0he03UaZm3QRg1yPFDoa4btM441xhhOzZb3Vbe97ecLIu8bNYin7h4FnW9uMtlSZI5GOloujwVwyahjldi0Rmj159Ew00evEVXathMR2c4Vk10qu6O0LgAMxkaSFUQ4y7BoRIAp94wCOxzmurD+yTU7a6+iarf6XPMsrIIul9T0/X5ZXQgeD4iXUViA5JBdL12Qg5UqQxls3GeIMvBU/J8DtJxRb90Gdagwko1rW6gxTLS6zdiCKcTGZCtLpwmeLayQpvlWTqBirMk4+GuSKAiUUzn8r4d29nJeWMUVqo74q9DYI3cJt5N3PYzkjnhnU3Uos8e3isiqxqiQQ1ZG0qPjnGqNJuYNGGQKmg/fkEyqcWuJ9ytdwjhymPtl++puviWs1iPn6ZtP3c0Oz2OxzL1RAvqsf2Lcwp+IrVjKqkL5KTBhCY5tberFQOkwbPSOjiltLM+IpsJzNj2o07xFq/UcGTlgWgZOAY5Hbo2vDWT7gzTVZrWzbrdhaFslpr9TcSiCTO0W+gY6mCBYYxYkGkZZym3zR30l2PEtZSTtBaIRI0aDu3naRXBxjDEMCDkIp8o0tU6f07SgLXVOmtc0iNJWjOpKJr+9lby+YWJSG3cAbiscMgC8K104DSvP2nKNVe6agGL0dYG47CTt8wCA8CA8hx9e/trsB8mIFMUDHKcvWQSAA9Ze3cvJgEeQEBDkA/Pv21F16vOe2pdE2JZOT3VYCXk1Xtm/FWTH1tzrTFWpjkkmWPHbhjItMmOZVwumo3r07YaO2qyDAzdq9flfKCjtNhbc3jDO5pZjWnzxnc66Zu3tuP7KxLCXOlPlgWB5DTDAF3bY8nGuEBayJY5++ZoOCAUHp+tMT9C3vbe5XcrGJM7d0g2nd6ghsLuHfCswGcHBBx5fVek9UsLRtVsoRqmjIqs15DIRPboSFzqFmEaeylVsqySKUJwY5XR42bZ/wBekZRRInJjpdAnDjsHX1CXuIgHcCiPYR4/bruFwUBABAREft/3/wC+vPlEAHqMdwYEy+SADwBDEMIdJhET9zl6R79PsYe4dx1Wj5vUXyzETS57Ab9ce/P055Hj7D+Y+2tiSS2QhTdLkg5IQsueAOVDZJ5wOM888CvHo5lkMcc0GxFWR7hiyqSx5hjUqN+wcGQlSSeVBBq7CsTkQ57h9O3+v/trkFQEOeB/l/P7as51Tip0kTAQ4EeTAHIgH1+v7/f9vOqgDEEvQJQER+xu3I9+w/mH7/fWMPFx/SSc7Tt8LkDjOf8AtnPBzzkVZbiF2baLgKjbCxhO1jwPI2eec44wScZGM1XeaHIB0+/tyIf+/f6du+uRWKHbvz+737/n+Xt76tpEgA/UJRAPyN3/AC+nv/XPYNdTggqHDgSCHIfrCAG4Hn34545D7j9+OwayAwMwH0kgYBz4ePRTzgc8/l39qlmcqGhCk7hxd7rZSOOzKsufjyjOf33bzB4/V/Z7ccf/ADfs08w3/L/l/wDV+zVKBE+A9+wBz8wdv4iA/wAQDXPQl+Y8+w9Re/7O/wBu+niW/wDf+3oPj9x79/f4pun/ALqH0/8AcMfbt9lyOe/5/Fd5lxKIAJR7/kH/ANWuRXAB46RH8+Q+/H351SnSTEOeoS9PcB6g/wBeP4/y1TCCZxE4q8Kc8gIGAOoA9x4EQ7ffg3bsHHfV18JhvE42A45Ud+Bydw455Pp2AqU8WRijxtAO6zRfbofdXBClcehAPyRV3BQohz39uf676plnhUf1imN/+QoDx25+pg9/YPzHv276sElPxEO0VeSUtHRzRqTzHT2RdINGbZPqAvW4cuFCIol6hKXrUOUvUIB1cj303z34imybbJIwSGe9yuOseHs6DpevpyUi5k2koWPMim79M6gW0qg3OgZykVRu7M3VX6/MRIoRJQxMEt1awffk3H0Vc5+DjHA9ySOPStiO2vpJVgttNvryRsBCkLKjg48wYBi2M5wFyfat4E5RJQS8JLF5J1HKchSmSHjkE1AA4gU/04L1ByA/Nx318qy7RIi5l1PSkRTA51XAkIQgCQTGEw9Zi/ouP0oGEoF+4h31DFZ/F8w/c4zJDXaJRchbpbRUo6mvK9YaWxr0Riaad3chnUVGktcvZmk0gSMj0ni807JVFztfTCBEVzGDUeucXHid5KnbAz3S7zdsmyXblN2pw7pdl27P5TNeSJGWUM9+DYwtkN8Dx0M3WpmAVl3Nhl3ksVMzmGZtVI1cJAVG/Nm1eAGMp4iIwbL+GWAA+84LFUZVwRneOeB6491pf7O9UkG/WJ30kOYzHaXdrKdWlEgQotjpsIkmuS6NvVpTArAhtwyKmH3c74lMaUq/1DbRUA3MbrI+tHlKJiGruBM2Wcv/ACyMrJZLImi5YxFZj03BTu/h5Jh6Lldkw9ECTlZ020922+LflDHu2yvSfiH7b8h47zjTjJweU3uKoup2Sij5DlvGN7anHntsPYoUZJYwupesNICTNXFFiskHsokkZyGqOMa/tUp1PbU25+JbkeKrVcfuk6RDbN6fZ9i9Si1pg4uZ/wDEldxrL21pcJ506atlIyce+iXhmyck0bkXJKqmTzvQI3whcX2+v5kwxtth9xe5akO3N2eZdpuKmOQtzjt+1ZukrPmC6ZEt6dSlZycdKPDqWu1Kyy009k5gFPSuvVuDpab6gkgX6NrkbnsALOMk8AjMYlZyc9+VJHZBzn0Y6Qa1SFf/AA76hvYS6st9qN5caTJexkrufCobOwi2nI3tcmN+HmlHB2b3H+LjgCHpZ61gh5O5uzVcXpavT8e0GOcozqDmQZPVHFlknk6MK1iYSFbIKiaTIs6ME0vDR5kEyPzO20duyDxGt2eLrpYMN3vwrNzDvCUIjaZSNzzQHFLtua8u3hxNR4GvWTcUrWiHha5M5BRVk7HcXbTIlnPHT/kMiqyhHRnyWZ9gee6J4vmZN22RchYmMTG2FjNsFYzrlujHhoG74oyZIK2l/K2qIloxq3kZ6VeY6rUk5jul61qzkFY9pJyaa/qtSWB4a+1iGKZbFdJNgGfE6xFrjgF9/ZfcX8e5MC0hFPrFX2yD48XIO0mjx80MZRNy7ZNFlAE6BBC8DX74mF3BcxOADJDbSwkdvKY3kkywO1gWlAZQvlUkgzrEfQumwT9ONNqujTypDNfRwtb9QPFeGIOLea7i+qSyW0crRMqWUu24eYM0ixIaxw18SOaJ1Jl8OPxNHaxBUV6f7CcXKqJFXHrKn0mz0XoIAcgmUpjgUoCA8dg1c0fGL8PqHSNH5Dzm3xNeI8/orjjHItZtEdfMe2NLkkpT7iwhYmfiWdjgnZVo+WbRc1LMUHjdUjeQdJgVU/oibTc54/WKGB93GSYNouAmsRc5xSu4148V92qdXmrVZ6+5qrNuBlSyTJm2cJySotFlTkFomBvMPnHiYUR6aFhqzt83DMW6xFSZGtdwmsHyL1i4Koq7ij0Wu0PIscxCMclaoNZFOxOFZtIp3rltHqEK3P1VujGN1wVaMnCkBlk3NghSieJJz6nYEHqwyCfDwdKWd3GINI6s0I+GgZItWuJNLvjHEeVCXca2RznCot885ypETASbPd0PxYPDzyfb69QqHuioU/bLW9GNr8QDO1xoyT4qSqxkCv5mvx0W1EqSKhxUfPWqXy9IKCYxSjuUXL+KzlA5Mj0U4e5um310RTD25P8A+J8FDngvPI/MIB9edR0ZA3EMci43suFt1uxPN18rdljgr2XoFHH9VyThKZIzWRcTBm0nPWmDeWiiIP2ZH0TKSFVjnb6PbovFYVm4D0oaZ0nbV/s9l1nYuo1bazsqfTllerxrRingGDi0X3QVV38PUev661aJEKVqKoLuFUyKeSHJQOYpRr9YICAbeYDcFDEEBskLldwUtyfT35GMCol6L6hjjMi6bd3ECrva8sBDfW6IED75Gt5XEasoZgZNu5Q2M4OJ/oW81KyJLrV2ww08i1MUrxWFlI+VIzMYomKDszBy4BAR6RApTiBjj3IBigIh6FN6gqQqhTCBTgIkE3AdRQHgDh3/AFTe4c8G49wAe2ohYzwy9kMcqoy2t2G6bOk3SRXNkr+xfMC23yKvJkhKVnJ5BiMeFTJYHtdSOeOg5V8msdgzeOGqIJkdmDXpg8MKsKNwMjvf8TFdNThVVIN7WTXAnEwifrRUUOmYhuvgwqgXqEoCAl+btumSB/6uUA/4uMdvf+ea83JbS2uTdxXaADJ/o54HAHBbOeCf3nGM1KiMg2KYSGOJTgUxugQ+boIcExP25Dp6jAACIgI888a5O+RIA89QDwAhz0gAlEQDq/WD5eRAOfuIffUVwbWfERr5zxNE8S+PjqTEpAxrDLIez+p5MuzeCalBCGY3PJEtmFhLXywJMSIoTdyko1k/m5AFpJdgko6OmS1SeGfFLpMUrcYXd9hjcTYYcgLM8D3PbJXsEVXIKoKkSGClM01u55Kn6VGNgMaRRnGGPrC4XXj2rE0akm8VctrLCxHE6HBA4H5Y/EMH37+o9xWsJo2AeJHmU8LGNyzlwPN9mVAEY5y7MoGc8CpaU3HmD+oYAD3EQAP/AOQ6ol3ySLkqShT8GTOcFQIUCFEvH6MTdQGBQ4CIlAC9IlKYROAgAGiFDNPjekESDsZ2EgY3YRDfXlU4Bz7CPO0gB/Z7B+7WFr5vd8Q9lFWPHt22Zt5HJidvrVdsUpsivTrcwbGuPJyCmJOcn50mQKxgVAlxVeMYqJga60duECt5N9OO5lg+hG0VI47gTQqGjljYlgCNhI24BJ4PBz2JOOeQBnHX0vTJNSlRZDbWqbd8xlnRlhi3IpmL9htLAEDLEkIoYsBUjuVMv269XJbBmAlkhn0DNRyvktcgrQeIq/Is1V0UmhEhEs9kSY6kCQde8+PZtY74xNupxu+hW0RJ4mu1wgdr0BWdse2SrI3rcdkMH9nh4WReHWbx55R8irec951tCLVw6axBZWSTcyckWMkZe0XScg2RY1BhKP5qKigtnilZtwNn3absM2p+HNmWoSO5hzc0iZG3muAxWeHtSJxfzWQHrmjtcxoXeQdJhJzc+9sD+uuZmxmbMlTmTk3Ei0nf25bd2GDIuxSr+xPMh5iyXKsrHmfNM6l6ew5Im45Fy1iRFmK788JT641fvIygUMkrIRtCriqNeiXq7RoRRSohW5UbywxtLbSVDEbTluAQcjB459eNoG1q2qQ2QXS7BZDbiRWB2hZ9QkTKpLMFd1ScKx2W6s8cSswDMxeR6nb7tmZYPj7NNy1mkMoZkyW+ZTOYcw2NsVnNZAmWCK6UY3JGg5kRgKbVUXTmJoNNTlZJpT66ZCFav3aTUqx9l0WigkIVYhSB7mTIInKBjcc/OJS9QB7AIlAeO/Aew3Qo8gA651siVhnGBn2Hb2xzxivPoMStPj7Rx5SQPsQwG4RceTcQC3fJA7YqlRSMQ32KUOCh/D8g5/b+37hqq001QnJJPqc/voq7c5ZmLMWZmOWJOMkn9P5700001FWpppppSmmmmlKaaaaUpppppSmmmmlKaaaaUprg3sP7B1zrgQ5AQ++lMZ4PrxVOUOBNyPID7e/PPf8An7c/s18HEAAeAEDce/AfyHkPp/1DVR5f5/1/EeNfXQHHHJuPtz/X7dHG/gMyg7cleDx6euf59O9EQRIUXc4PfLbDk45yM9sfPOe/rbyFUE49RxMXgeAH6fYfcf2e3PA67ugf/T/Ef67f+/OuXA+QTrKY/IGDgpS9ZlPf9GHUYOBMI/rc9vtq1GkF0+kyhyiYqSh1kEQ8zgw9PlkIoIlE3HBufkATCI+3T3xoqqdhndm4PmJB5I4HGPQnv3qkcQLLGkKksSTvIlK/dJJc4wDzgYJH64q6dA/+n+I/12/9+dfXSX2EA549uRDn+vuOrSMgsBAMBwP+hMuboAOkhOS8FEwiUeQAREphKHWACAgAl4HFeTs+Ynw7VZ685JyHVqjWKtEnnLBKTkm1aFjYxMSgZ6u2807zyFBPyQybc4j0j0lN34xXICx5eF5lDLxIAyAnG0kEjGecH24NbFjbz3s0kFnbySSLItvutYS26ZiqiNWQMCxbIwcZPAJzzl52Q4kAC8hyoXuXtwA888n5AUw+gnABEOQDp+bnXmJOSj4kir6SdR0THIIHeuJOQdpNm7FVMSgZ04K48pqmmTrEDOVHJOkTB2DntHKz8Seo53m7Rj7Z9SrJuBtMQwr6yVxZNDwuC453bE5FavFt2RFPVSkZHOEoqQO6fQ1QsPpvJTAqSwKgIY+ybsJzLvCOzbbw9wNoSxVY6NZ61dNq+HJGQq+P5IbQ4g3KsdcLy2dNXOV4itfCjNYh9P02IWXTkHiosmgrKJjoNLaJgGyt9/OzbGGCtkLksPKQrcYUsQc8V65OkZNKO7qTUodHkkVXWzvJ9186soYBLONZGRypyBdG2QnCbwxOKbcl4wOEsZubjQcCIjuLzNVbI8x7Js67JrNsR42yG4OqFdhs731m1lJ2nsrL8PmPh03UqTkFEhoh6CxEgOiK2ly1C8QffEWKms742m3lVfAK142pZLsrvBWAI2Ps5QVlMdT6VVY5Ib70arHHYoh8XyVUsVKswbpAlDHGwuwj5kNtnh07Q9qTesKYVwzTqtO1Wop0iOuKUNHnua1bKDYFY6QshG6T56m5MzbKORUEAWURTOcvJQ1ugWKZk6B8sDGTSFEpzgBlAIbp6w6x7/OJSib7iUoj3DnWz4V++DLehwyLgW8YhCZHmBLqS/fABQADOSxO6sj3/R+knbpuiL1FIGybjXo3trWRhjGNMtbp9gU5IkF4JGwmVXa4ki3xpsHeIUem4lzBkhCexNRK20q9MwnjGkMsK4ZaU5BJJJSpz9BgJqbhLfENE0GSEWhIIMUotIrlNukcHqoJbq4gwbirA9VLRsQY7q+NqgnIyD8kBV4pnAxgSL46PrX7aNYJFbg4kRSTO8WESnOKKXX19unPJWKReRAT9RgKUxxMIqGAv6vUf3EQ5HgR9uR419ejQDnpIUAMJjGAChwYx+OsxuOAExhABEfcR1mFpDtCSbpu2WmYyuoB+6JHy+wY4XcRzjFczUOq9b1aAWFzKllpasWj0/TVjs7KKRgA0kFrbpFBEzAckIWPOScktamyPJ/nBMwAH+EvAgH0+gBz+zt+4dVJE0xEQOUvPIgAgUAH8vmH7dxD948fa5lTKX6B/Dj+Qa+TpFNz3MHP0Ae38P8ApzxqWghjbfBbDeSBlXERVfxEEA5/+05zj0rzioyDw0nlKE8tcM1xLyV4DlkIHB4wf1NWh21IZMUvmMCg9AiYhVgIUQHuJTGKAFD/AJh7hz7a0iu/h2bM8h22eu9z254nnbRY3K8jO2KTpkI/k5109V86QXfOFmxlHB1FypKHOoocVDdzAAhyO+AtUxAQER4EBAf2CHH311JsEEk00ydQFSTBNMeoeopA47APPIAPAc8e/HfWObT4bth9MkFxEMFYZLdCFbjnfuOSO2doJGcmuvpeuav08Wl0LUb2ynddrvbXU1uSOOCyNuK98oSV9vQ1C1kvwko1lQ5qp7Nd0ecNj7mayXI35c+L3jyboDU0+d+7tTJjjFKxVCGjiSkipGuGr9GRUVi0WR2iTdwR2dRGN3KWKt/2HZiLbbutwFphKLjRy8iMT758D4kaWexxuII101Tsln3TWha4Viz4xb3F22qDyVq1RTyq3sb9BRZ++RNBtDvP1kizTHrAx1Dgc5j8HOJgL1f4C8j2TL/hL7F4DjRNkkmmCYCcShyI8m5E5h45OcfcxxEORMI8iOsD6Tbyja0kiJ2McbOiY4OB4boQAeRhhggHGQK9Ppf7SOrNKuPrL6fb6nqYAVZdV02x1CJEKBHAhvIbiJ96ZSXxUk8ZCwc5Ymvzp1PxAN3YY8f36nyW2HdRiBvkWwwyGZILJtsxvPRdNj13BULDdsWExfMyEXDQQFaM5eUZv5YTvZVkCCTghjKF3Yxf4j42un1t64237gLrLu4yMdyk5iGuV68Y69c7bCqZ5XbLPW6nWCTrjsxVFIaSmKlASLtsBVF4dmoJ0SZw3B+H3t8zVdnObGkVPYh3JOmEbBDuRwlMK45zY/rUUmp6OjT99hSJTk5jtZyVk6laY5eEi5VWLjQX6PTJmLFbk/aLv620XCUzbjyxwO5BSWkpC7ZYlcOUtntvyNJyzFb/AIRM23utSVspm6masCkk+O3G/wCTsbhVEEnyjJeSPKrpo6P1fe2TJDZ30iWWQ5DQxtKsjHDbXwzOucD7UuwALNJgE16WLqnovq23mtep+ndC0S/WVrh9Q0lbjSNLmkkaMgTWUUzxWgba7gWCpEhbZFZKMA7aZH8ZXaLjPJuKMX34mVYC1ZPkpGHdt3dMZJI4mlWHkefCZUEthM5jXhzqnSIWtNLUkKjJ0B1UwIkZbNw+Ift5mDGY4sWvec7GVRUVaZiGqrSVoj49E4FcS7lhYXlXjjxsSqduyePCyhnRF5BsVFosRVU6MU9G8Ufw7ci5qxeXxAsT44wfvPwSxs1nqcpkWBTsrvb9IncRLax01hkGXiIyShcksHJmbS2Q1bj38U1esOhvYZAiaSikiDHxaNmd2at2e3qTu+8a6qenkksX7a6Kpb8jpVtMhyPbY5Z3OQoVfGuQa7iOZP5FCxrvAcy8eDePcJqrKIZYrzUTK0AvVYKDh/AKnCgerSEKxJByyDudq+g0dWttG0tbeX/Y/WzprIBBqEOvxXWmXRYswaKSDTEBiVBjy3KuoClnLMWbIye+ldRY3mbPd6xBOJSotk8W0rpbDwYA884ZRDzDLByJC9IAkCZw6jcjq3y267Pd0RThcD7ScnxNsB8l6l5uXMxw9j48G2TXTdrxttqDrKco4sYuFGXwmIVrbZq/ZHkl15VkoyRQd4qldyviV5PerT23LYnjam4ncK/DU3W8LNk7g/OEbJoiLacl3+KKPi7MladQLZRZN3WlByL5lkSTUI+bw3ACPoy7SN9dq9PXcl+JUpK40lTlbWiLw9tYrWCMhvobpMLgtSyzXMuTU3Q5IVioGRmo6IfOUUwOkVuJVTCHRjXVGGFvkKkY86Rk/h7eUAH2OOD2FedfX+nVZWtenWkwwwrX91fw7xyVdSUIIbnbv7HmRWGTQTeafEhrpGv4jwvs5hjSIqKNwmd2NoiVHaTYSlMmQJHCDcz9JsCxSrnL2A6iXWUonDjSCx+LDZ071YsSX7cFtP28WKoTMhHO39GsWRtw9jeWSGcGjS0NpTpXEtAq6Z5Z04MBbCyujxNk5YtmiCLpvJKvGu+n+6Z2X21cA3IVq3b4V2BCGqcpvgu73crK47ZrgB5dvj5a/MnBqmwsCyccvPt45YhJleKilnJeuPR42hm8obZ9oOMatVY93RsbUCjtorHtLpFWbNWrGvN2zRVCGqUJXoJNweNTIjHi1Zs/TN2TbySIqqocp868kZj4a41GVe2I76NB2GfMtujr29HUgDg47dGy6nuNQuFtNP6M6aN4zgRhNMvrlg+V8jR3F9PaO3PKSxPGxPKHIqGB3nPxjcngStbd6zkyAv6Qndrvd5O1uhbfsQnhWX6GRSbXug5azfNnt5XDll8JrhqanHyrP4u7czserGNm7/0+RsPeJpc9tQPd5HiKY52YXWHsL6dkH+2eijZvLq0X1sGDRrZlZzHtglW8t8TRdzMejAmTYumLDoM9KbzEt3kMz7pNxlmnq5iI1Z22VFqkodGyZHbw17zeSbiFStZmHmMCJyKELFQL9Z2RZpagyG4eiLIqa8IgZz8mXMdbOKWwssLk3LT6Xzfl2vSKslBXq/FK8QpE09IcJ1bEtfdOZJLFdXl1BTMFegX7lqRJswQFQ5WaZtan0dZpR/SbpVGVIlmkuGDZUthmYSbgO2fJk8H1X1MurNpURbWB0pZSRAPJDouh6PPq0bAFRbreQwLa2L7sBzGbieMZd4SCEaHKb2zT25DCUdhWBl95m4xvZ61CwmQ8k7h8y2PGu3fL9ebNkiyd4rdejV8jqyQz8kVhZa9TJWIj2j1kRYXcuwWZpkW2d2VeBztJ2tUWRq85Gr5e/FUgnOWKOnWLeNxxOSrMHCFcmnuNfPlIhK11iMkZaKhZ0Xi7pCNl5pumKCb1VM85DaDj26QooNypJHBQVSkAAKsoqYDnWUAAADrCcoG8ww8gbkee+q9KLZAn5YpFOAFOU4qgBzq+aYp1TKmHgTmVOUp1BH9c5QMPfWxFosMUqzvJJKVQqryN4koD4yM7Qu0+zhyB5VbDOG4mpftb6iutPfRen7m70TSvpcN29yJBNqk9xANsTG9KpLAgBbfHa+BC7bWdHMdv4Meh/DA2DJCVu22rYO6Eh4Dpx5X0yIcAIlHoKzEh+lMDlBQDcgA9IAIGHjK2JdlW1fBFm/GuGcG47oNldRKkCawU2qxUG+PByCrZ6uzXWZt0VDslHEcyOZMT9InSTHpASgIbgA0RAhCAUSkTDghSj0lKHsHSH04DsHH0EQ9h1yLZPpEvfoECh0CPyAUnsBSBwBQ+/HvwH2DjcNhbHBbLcgnyIDkEc7lUNk8557Z75xXl77rjq6+hNtda/rF1bsCskFxqdzJHMpxkOjuyBeMFMYIPfJNWgzE5ki9TdoZUBHkTkKYAKA9uPl+39BqsK2SKQpTkSKYA9iAAFH93ACPbjn6B9frqvOmBx5ERD8v6/r311KtiKiUTCb5fb/uP8NbapCNoEcMY7krCMjtx37cf/Bryo8SW4kkmluEjdQFjW4doUIHeOEY2EngkMao/KTIY4AX3AOopQAA4+4D9+4fT6/lxrvDhMgiACAD9uBH9/t9Pz7fX2DVR5Be3f2AAAfr2/f8A1+4NffR24ARAO/Pv35+/cNXLgYCqNvBI7ew7AYzx3we/JOKxhWZizqhlXCpcsoaQouNoKgjsMg85P61QgZMTe6n8Onv+XBh/r+fyoYAUKToMI8c9Qj+z3Hvzx9fy+n3rwT49jD/Mf8x18iiAm6uo3P2+n09w+vtrEqRbyzIDwRtPK54+Oe36+tXJlAAcLdduCBCFOB58DdkjGAPX3q0uyp8EIoU5uo4EMco8CQBAw9ZjcgIB2APz6g+msd5Bx1R8o1eZo2QazD3CoWFsWPmaxYWCEhDyqJF01yEctnBFUVUgOgVTpOkIGOUDDwIAOssLNk3BDJKgB0zAAGIYAEpgAQEAMHsPAgAh9uNdZmSZjdQGUL+r0lKbgpekBAOkPYOeREePr+zVwkBBVo1Csu1uOOcA4AHp8555781nhvb6ycS2zhmVlaNUZrd4iMYZJgzncp7EKpyMgjio8F/Db2oQi5JDEtGW272cA8ha6benwYjt8nBD2PW5Ow1lum+eVldcGjxxEKj6dZ/Hxrk6YHapiWmdbVc648KJcJbvcjV+LdD1zBM4xSu4h8u8SOVNiWBlbNZ4JxV2HknVI7YM0XST5YW7pU5DNSENI0Vk3KIG8sBOBgP18B1icCiXzBNxz19JjB1e/Aj30BokAAUhlEwATCIJnEvUY/PUY3H6xhEeoRHv1cDrQTTreMna8hBPZiR7ZI2MpGcc8855B5FejXrTqLINzc/WWPw6sItSQ8YAZLuOVHC/hDKQvdcVGU+kPEuoT5SAgKpt9z5DIkKupkC05FsOGp94oskZZeNNRIDHV/jGyUQ4EI9vIFsx1ZJqUz1Vq0UMLbV0Z73rdWW7eLyTtI3LsrHGNkyXKbpNOrloxdHvWyYBOva3b31xgZqdqLZwVd2xmnVTi372LRI6Uhmy6gtCSOEjkCEEoGVOJ+nzFFDidVXoL0k8xQeBP0l5AOfYB49tPhrPrA/kkAwABeQKAclKAAUo9uRKAB7CI9+NSbRWOGZjGAAqptjIwfxNtZ3+d7FjjliTmts9T6bejdq/TGkyM39Yul/SNKMzYUGSZYppbf0JVLe3giTJVYguBUG+KPFhx/vrzZkXaPs+PcoPJVJcPnGRsmXiqt4eAhcbMXgRTy3YyVYy8spbpuQevIhnX4ueCluiRco6sRlm7yHJFOpdcQYpqeI6ayqdQauEGqDp7IychJLi+nLDPyi5nk7Z7NKnKVxM2awSJ1ZOemXXLqUkV13bgxlVTDqgxhtpwDhWbyRZcSYhoOOp/MFvPf8AKMxUK3Gwj++XZU79RS1WhwyQSUl5xRSUkTqSDwyjg5nrkwn5VPzmkrcheQAR7gID9x59+R/aHI/fv99ZLaAxMS5GCMYUHtkcMCzbhjuSck/AAHH1HVbWWNrXSdNTS7MkyeGJFlklkxgfSZljh8bgkKyxRqgyRGGdiaU6ggIkNybgOvqN37CA8AHv37/f78e3enQMdQ36UBIUph4AvHSPI9u3bkffv9RD9w3H0xBAoCJhAvPACIiHceeRAR4H94fs19igUenkR+X7ft/brf3oqkIAGI+9gjGPb5Iz8CuDEGhhJA8S4Y4UyEEWyHGURzy3H3mwCxAOBxjtKHAB9f8Av31zp7aawVcZwM98c/nTTTTSppppppSmmmmlKaaaaUpppppSmmmmlKaaaaUpppppSmmmmlKa6FHBEzdJin7CXk3y9IAbnkwiJgHpJwHWPHbqLwA8679YIz/nXC23Wh2DKGcr/VsfUuqwclYZmXscm3ZeVBQ6aass5QZCYz6VTbgo2Ks1YtnCnUsiAkDrAdKVdcx5qxvhepTNsyLboepwsNEPZt4/lXaDfyWEeQp3DlFA6hXLoUwUIAkapKmKJygcCdZeqFPFHju4Hy9lWfotJxPlSeQRhGL6pqxkLXEpV7Ixqj1Kyqz7eTtMewioCIO5hDNpBvJyLl0i/dnGPbi1IVeM6s+IhtN8d7I7F7j3ZNnPMOdtqd5yFE4IZ2ifqda22q03IytfLXLTuSl07BOW+vwk6elKGscbUsZZLSq7dokZutLmkRTSkjx34H2JMiS1WyVu/pWJwcRblrOw+2fBNUaV3CtJF0Ins+PchWoWcS+3P1dXyItKKkbpQ6MoQreQMaHSGRVKTzt8OqG1GNNPksF0gqomjkRzeu7DEnhSDyRqO67gTu5J28V9n6I1H9iNh0pq8PWnT3UvUvXN8sqaG9ndQabo3T5jZTBdysPEutSlnIInjBtI4YQUjZ5HZ68ncN82/HO9guuJdtMFjgmYKzIuKjMNsUNXebcP1q0vlDlYtM+ZctUNjeUxJCAmxkSQ81j+mZc9eckiuq0ZjHNQf+mxn4I1XzRPwWZ/EPsp87XtJ43sDDEiku7stSqyLwDKyOLsh5MmUG8vuWokMoDZvVV7hTKb8LQJJGQh0vi65EZwcc4xpWKKRVsdYyqsNRKFTopjAUylV1gjDV6rVqIR8iMioiJZJAzj0GKJgRSbNylSKTgpeAANZLbgQ6XWUgF7m/wgX249g4H39wERD69g11YNPjKl5VuhMwAk3yEwHDggmPOzf5TtKqDwQeN2fBXPWd2Fks9Es36Zs7djbyLoSLbfWKyLtcX8uXuJi/mYxy3E20uwEgj2KuPsfYxpGH6PV8ZY1q0DR8d06GZ16q02sx7eJr1cgY1IEI6FhYhqmm0jopmjyk2ZIFIiiQAKUvGvbN2gtyEKgcQbABCemL8qaZS89IJcDwmJQDg3Bfm7ciHAc3E4icvSUoGEB4HqL2D27gHfkOftx+X5/aaZwD/CUwd/l7B/kH5ff9+t4RW8SrlVDZAJOMnGADwAe5POfevKGW9Bbelu1s+S0B3LO24g5fLMN3fnjsAMd6qy/qh+zj+HbX1oHsHPv9dNRQdhxjgce3x+lNNNNKmmmmmlKaaaaUpppppSqCQbmcI9BBIUwGAxDnL1+WcAMHWQnIcnABECjyHHIjz24Gwt2Lg4iRRYweSqZMqTgxnfqmZRECGcKKdA+oVDoOc/BxREolKJwUEQ9Qpz244EO3PP7/6+2qIROZXgolEoD2A4CA8AH379vf8Abz3/ADyoVKsh25PHfBIOOD+fpnH596154432NLbNMEbMbxYLxscZMqkjcgI4GePbJzWj+7Dw69pm9yVxhObnMR1jJ0vhqzDZ6NIz8aykFWhFklE5SsuQdt1hcUqwK+hf2CqicrCVkYSCeOTGViW46x272y5VwaZw92kZLXb16PbC7bbcsmrOXmO5JtFmKhC1Cg2kp362EKnHNXLlNdrXKnYE3QJxpBaFKyTNqSg6XUp18Bz9eR+v7Onj+eqBdmUy3WZZUSHDpOkoHmJKCHIlAoCYAIPbkeO5uA59tYXt18NlhGx32kkAYIUcq2O68k4LbTzlSCQe7pfUur6Yoto7+2j0qRwbjSJ7YajZ3hyPO8cjItrIF37JIwZYmYvFLG4BGlmLN2VWnJkMX5acxWLc9sYheTtmOJGScrQjdrGKtmy89SrG/YxrOepcg4eJmgZV4hCS8wzFN25gI86Z0CVWQN4mJqXPL1SqtbVl/I6bZtKvseYcgWVps6EXKFMqSdUcPJaEhxjWokIi+WJKneEWeNiiyEDnMT63l7KsQ77sTJ4Wz4tdwpDS5Ql4bfgKzvKfMlmYRlMsGJVZhkVRc0YKE07M5jwTFMVitTgYwogOsk7c9tOLtqOJq1hjDkW7hqfVo4rdkd6t8TnpZ8YiKclYrNKqFRVnbJLqJJuJSXdgV09cCouqJjCOtGSzmjORuVSoHlckFiwBXGOBgZ3Akc42jJFd6PUOiihuY4LuXUmfMmi2cwg01125DNqiqZ2G9tgsTagqik/TWJGdZWEJuz3MOX7i2zqW2LCzt0U0dWq2ZeSzZkGnvhVWbL2SbKWDNg22Q6STVI8fWX9/brKSD1L4mUjFM7rNeJNnW3fCNqUyNSsdxH9qEmxdw0zl6VaN5fJ1tCZVbvZh/arguglLTb+aeRzZ5MOnapjvXiaayxjnKAhtL5bgrgfL6ATABKfjqO6XSQ+UomOYpQE3BuwiYQ5N9eR4+hMCyKhymdFMoIGPwAouUCKCIpkDpOJQMQAEoiVT79uNZ47OfAWSRdhxlUJwB6Z7ANjjKgbuSQxJJ077qfV5reWy02SPTtHlDJLo2jsLWGRCAGi1KfBmvp2OGzOTg4VPDRY40/Dh4nvg4b8sa+KvRd1/hLWLMFfsu7OXu0nnq3sLbJUWo4vlX1jrdhn4+3ZHhXMnPt6Hkl4VZZWEbVZ2hChCtkE0pIqnmI/qowTvZoFrspcL59Ur+CN3NXg3EhkjDMlMuVodu3jl2bJe0Y/vMnGwjC6YwmXjoqlKsUi1rVinIwCPX1PhlPNbpbwmBUpllTEOYwkDrFsoYXAiHINw8oxSJmOYonOoIrABQLxybq1hLPe2bA26Kls8d7jcP49znSGU8wtsdV8k1iNsldJZo9nIMG8wSMk0HaDSSQZysgg2dJpnUKi7cp9RAVNztLbwxjaoGQPYE8EYwSP0zx2x815NRYYSV/HtkQKptXD7lIC7sjcQxzg5JywHGOazRCXOqWYjg9asMLYk2aqKL08FKx0qRidymqoh6s7JysRAFSIqGIBzAc5SCJSiBREPStV0XKYqoKFVIBzJiYvcvWmPSYAH68D25DkB+giHcdVtv2zzaptLjLXE7adv+K8CRFyexb+5s8XU2JpTWwOolB41hnkqhCtEiyLuLayD5qzVXATIIPnRCCUqpwHaKMaosmpGjdMEUEBEiKJREU0Ug7ERSASl6EkwDpImAdJCgBQEQ0KMuecAn7o7AEDBz/pVlk8Qhoh9gQQD27Y28cDt8HmrhpppqtXpppppSmmmmlKaaaaUpppppSmmmmlKaaaaUpppppSmmmmlKaaaaUpppppSmmmmlKaaaaUpppppSmmmmlKaaaaUpppppSmmmmlKaaaaUpqBvxBPAlwp4le7ms543WZvy/a8NVTGKlNq22CBknlZqlbsh1k1Ze+wVuYzpnsZI2MiLBKfjmtfTSlAiowzt4qLJACNNKVm3w7/AAUtgXhu2+z5W20YzmGWQbrEpQq1uvdmeXqxQUB1LGcw1cmJlMzyHaSxlEjTxWaiRZkWMYLwhhYICEuvSX7f56aaU9c+vv606S+3H9ft99cgUoBwAcB9tNNKVx0lD6fzH/XX1pppSmmmmlKaaaaUpppppSmmmmlKaaaaUp76+egv2/mP+ummlK+tNNNKjA9h+4VxwH2D+AadICHHHb3+3+WmmlAqjsAPyAFOkPsH/Xv+fvrjoL9v5j/rpppU1yAAHtpwH2D+GmmlQQD3AP5gU4DjjjsP7tAAA9tNNKAAdgB+XFc6aaaVNNNNNKU0000pTTTTSlNNNNKU0000pTTTTSlNNNNKU0000pTTTTSlNNNNKU0000pTTTTSlNNNNKU0000pX//Z)![ref3]![ref5]![ref2]

<a name="br6"></a> 

6

Fast Approximations and Coresets for (k, ℓ)-Median under Dynamic Time Warping

▶ Deﬁnition 6 ([[3](#br24)]). Let H be a class of {0, 1}-valued functions deﬁned on a set X, and F

a class of real-valued functions deﬁned on R × X. We say that H is a k-combination of

d

sign(F) if there is a function g : {−1, 1} → {0, 1} and functions f , . . . , f ∈ F so that for

k

1

k

all h ∈ H there is a parameter vector α ∈ R such that for all x in X,

d

h(x) = g(sign(f<sub>1</sub>(α, x)), . . . , sign(f<sub>k</sub>(α, x))).

The deﬁnition for the sign function we use is that sign(x) = 1 for R ∋ x ≥ 0 and sign(x) = −1

for x < 0. Observe that the class H of functions corresponds to a system of subsets of X.

▶ Theorem 7 (Theorem 8.3 [\[3](#br24)]). Let F be a class of maps from R<sup>s</sup> × X to R, so that for all

x ∈ X and f ∈ F, the function α → f(α, x) is a polynomial on R<sup>s</sup> of degree δ. Let H be a

κ-combination of sign(F). Then the VC dimension of H is less than 2s log (12δκ).

2

Theorem [7](#br6)[ ](#br6)implies a bound on the VC dimension of range spaces deﬁned by p-DTW for even

values of p, as follows (see Lemma [43](#br26)[ ](#br26)in Appendix [A).](#br26)[ ](#br26)The decision of whether p-DTW exceeds

a given threshold can be formulated as a |T<sub>m,ℓ</sub>|-combination of signs of polynomial functions;

each one realizing the cost of a traversal. This yields an upper bound of O(dℓ<sup>2</sup> log(mp)). The

situation becomes more intriguing in the general case, since for any odd p, the cost of each

traversal is no longer a polynomial. To overcome this, we investigate range spaces deﬁned by

approximate p-DTW balls and show that we can get bounds that do not depend on p.

The following lemma shows that one can determine (approximately) the p-DTW between

two sequences, based solely on the signs of certain polynomials, that are designed to provide

a sketchy view of all point-wise distances.

▶ Lemma 8. Let τ ∈ X<sup>d</sup> , σ ∈ X<sub>=</sub><sup>d</sup> <sub>m</sub>, r > 0 and ε ∈ (0, 1]. For each i ∈ [ℓ], j ∈ [m] and

=ℓ

z ∈ [⌊ε

−1

\+ 1⌋] deﬁne

f<sub>i,j,z</sub>(τ, r, σ) = ∥τ − σ ∥ − z · εr .

2

(

)<sup>2</sup>

i

j

There is an algorithm that, given as input the values of sign(f<sub>i,j,z</sub>(τ, r, σ)), for all i ∈ [ℓ], j ∈

[m] and z ∈ [⌊ε <sup>1</sup> + 1⌋], outputs a value in {0, 1} such that:

−

if dtw (τ, σ) ≤ r then it outputs 1,

/p

p

if dtw (τ, σ) > (1 + (m + ℓ)<sup>1</sup> ε)r then it outputs 0 and

if dtw (τ, σ) ∈ (r, (1 + (m + ℓ)<sup>1</sup><sub>/p</sub>ε)r] the output is either 0 or 1.

p

p

Proof. The algorithm ﬁrst iterates over all i, j. For each i, j, we assign a variable ϕ as

i,j

follows: if Z<sub>i,j</sub> := {z ∈ [⌊ε<sub>−</sub><sup>1</sup> + 1⌋] | sign(f (τ, r, σ)) = −1}̸= ∅, then ϕ<sub>i,j</sub> := min(Z<sub>i,j</sub>)εr,

i,j

i,j,z

otherwise ϕ<sub>i,j</sub> := ∞. After having computed all ϕ , we return a value as follows: if





1/p

X

Φ(τ, σ) := min 

(ϕ<sub>i,j</sub>)  ≤ (1 + (m + ℓ)<sup>1</sup><sub>/p</sub>ε)r,

p

T ∈T

ℓ,m

(

)

i,j ∈T

then the output is 1. Otherwise, the output is 0.

We now prove the correctness of the algorithm. For this let us ﬁrst observe that for all i ∈

[ℓ] and j ∈ [m] it holds that ∥τ −σ ∥ < ϕ . Further if ∥τ −σ ∥ ≤ r then ϕ −εr ≤ ∥τ −σ ∥.

i

j

i,j

i

j

i,j

i

j

This follows from the fact that Z coincides with the set {z ∈ [⌊ε <sup>1</sup> + 1⌋] | zεr ≥ ||τ − σ ||}.

−

i,j

i

j

For all i ∈ [ℓ] and j ∈ [m], it holds that ∥τ − σ ∥ < ϕ , which implies that Φ(τ, σ) ≥

i

j

i,j

dtw (τ, σ), as for any (ℓ, m)-traversal T we see that

p









1/p

1/p

X

X



(

)<sub>p</sub>



<sup>p</sup>

dtw (

τ, σ .

)

ϕ<sub>i,j</sub>

≥

∥τ − σ ∥

≥

i

j

p

(i,j)∈T

(i,j)∈T

![ref4]![ref4]![ref4]

<a name="br7"></a> 

Conradi, Kolbe, Psarros and Rohde

7

It remains to show that if dtw (τ, σ) ≤ r, then Φ(σ, τ) ≤ (1 + (m + ℓ)<sup>1</sup><sub>/p</sub>ε)r.

p

For this, let dtw (τ, σ) ≤ r and let T be an (ℓ, m)-traversal realizing dtw (τ, σ). In

∗

p

p

particular, ∀(i, j) ∈ T : ∥τ − σ ∥ ≤ r, so that ∀(i, j) ∈ T

ϕ

≤ ∥τ − σ ∥ + εr. We

∗

∗

i

j

i,j

i

j

conclude that









1/p

1/p

X

(i,j)∈T

X

X

Φ(σ, τ) ≤ 

(ϕ<sub>i,j</sub>)  ≤ 

|∥τ − σ ∥ + εr|<sup>p</sup>

p

i

j

(i,j)∈T

∗

∗









1/p

1/p

X

≤ 

∥τ − σ ∥

<sup>p</sup>

\+ 

(εr)<sub>p</sub>

i

j

(i,j)∈T

(i,j)∈T

∗

∗

)<sup>1/p</sup> · εr ≤ r + (m + ℓ) εr,

1/p

≤ r + (|T |

∗

where the inequalities hold by the Minkowski inequality and 1 ≤ |T<sub>∗</sub>| ≤ m + ℓ.

◀

The algorithm of Lemma [8](#br6)[ ](#br6)essentially deﬁnes a function that implements approximate

p-DTW balls membership, and satisﬁes the requirements set by Theorem [7.](#br6)[ ](#br6)However, it is

only deﬁned on curves in X<sub>d</sub> and X<sub>d</sub> . We extend the approach to all curves in X .

d

m

=ℓ

\=

m

▶ Lemma 9. Let ε ∈ (0, 1], and let m, ℓ ∈ N be given. There are injective functions

π : X<sup>d</sup> → R

<sup>(</sup><sub>d</sub><sup>+1)</sup><sub>ℓ</sub> and

× R × R

π

:

R

X

→ R

<sup>(</sup><sub>d</sub><sup>+1)</sup><sub>m</sub> and a class of functions mapping from

F

d

m

ꢀ<sup>ℓ</sup>

ℓ

ꢁ

m

ε

(d+1)ℓ

<sup>(</sup><sub>d</sub><sup>+1)</sup><sub>m</sub> to , such that for any

, the function

(

) is a

R

f ∈ F

α → f α, x

ε

polynomial function of degree 2. Furthermore, there is a function g: {−1, 1} → {0, 1} and

k

functions f , . . . , f ∈ F , with k = mℓ⌊ε <sup>1</sup> + 1⌋ + m + ℓ, such that for any τ ∈ X , r > 0

−

d

ℓ

1

k

ε

and σ ∈ X , it holds that

d

m

if dtw (σ, τ) ≤ r then

p

g(sign(f<sub>1</sub>(π (τ), r, π (σ))), . . . , sign(f (π (τ), r, π (σ)))) = 1,

ℓ

m

k

ℓ

m

if dtw<sub>p</sub>(σ, τ) > (1 + (m + ℓ)<sup>1</sup><sub>/p</sub>ε)r then

g(sign(f<sub>1</sub>(π (τ), r, π (σ))), . . . , sign(f (π (τ), r, π (σ)))) = 0.

ℓ

m

k

ℓ

m

Proof. We ﬁrst deﬁne σe = π (σ). σe consists of m points in R <sup>+1</sup>. The ﬁrst m points consist

d

′

m

of the points in σ together with a 1 in the (d + 1)th coordinate. The (m + 1) to mth points

′

consist of the point (−1, . . . , −1). The point τe = π (τ) is deﬁned similarly padding τ to a

ℓ

length of ℓ similar to σ. Let τ and σ denote the ﬁrst d coordinates of the ith and jth point

i

j

in τe and σe that is for j ≤ m the point σ is exactly the jth point of σ, and for j > m the

′

′

j

point σ is a vector consisting of only −1’s. Further let τ <sup>+1</sup> denote the (d + 1)th coordinate

d

j

i

of the ith (d + 1)-dimensional point in τe, similarly for σ

.

d1

i

ꢂ

ꢃ

2

The set F consists of all functions f (τe, r, σe) = ∥τ − σ ∥<sup>2</sup> −

r

, where

zε

\+ )<sup>1</sup>

/p

ε

i,j,z

i

j

(

m

ℓ

(m+ℓ)<sup>1</sup><sub>/p</sub> + 1 ]. It further contains the functions

) = +1

i ∈ [ℓ], j ∈ [m] and z ∈ [⌊

⌋

g τe, r, σe

(

τ

i

d

i

ε

and h (τe, r, σe) = σ <sup>+1</sup>. The function g has k = ℓm · ⌊ε <sup>1</sup> + 1⌋ + m + l arguments,

d

j

−

j

corresponding to the sign of a function f , g or h . Always ordered in the same

i,j,z

i

j

way. To compute g we ﬁrst use the sign of g and h to infer the values of m and ℓ ,

′

′

i

j

that is, the complexities of σ and τ, as sign(g ) = 1 if and only if i < ℓ , and similarly

′

i

sign(h ) = 1 if and only if j < m . It then invokes the algorithm of Lemma 9 with input

′

j

sign(f (τe, r, σe)), . . . , sign(f

′

′

(τe, r, σe)). The statement then directly follows from

1,1,1

ℓ ,m ,⌊ε−1<sub>+1⌋</sub>

Lemma [8.](#br6)

◀

![ref4]![ref4]![ref6]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEADADASIAAhEBAxEB/8QAFwABAQEBAAAAAAAAAAAAAAAAAAYFCv/EACsQAAADBQUIAwAAAAAAAAAAAAECBgADBAUJBwgREhYTFBUZIVmY1jFhcf/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDrftBuCDaStlItxvq1AEGKumkbOzpOza8dpZEp073KcZalZDo2O4PKCZsruC3uIylKUNqOGIxEqpkg+l0G9GoLU7IJ3JTZHd68SOy449CF0EOBfrEf1jGDQ5Yhe4RU+8sR9BZyxC9wip95Yj6CxjA5Yhe4RU+8sR9BajRdwQ1mC5Qy3d32agi908rpRMTpC0u8jqpDKAjgz0xpYpU/o2A4rKYnoWKg97cbUoAG0L8sYwf/2Q==)![ref2]

<a name="br8"></a> 

8

Fast Approximations and Coresets for (k, ℓ)-Median under Dynamic Time Warping

We use the previous lemmas to deﬁne a distance function d]tw between elements of X

d

m

p

and X , which we will use throughout the paper as an approximate function of dtw . To get

d

ℓ

p

an estimate of the VC dimension of the range space induced by balls under d]tw and decide

membership of points to these balls, the approximate distance will only take discrete values.

p

▶ Deﬁnition 10. Let ε ∈ (0, 1] and deﬁne the set of radii R = {(1 + ε)<sub>z</sub> | z ∈ Z}. Lemma [9](#br7)

ε

deﬁnes an approximation of dtw (σ, τ) for any σ ∈ X and X , by virtue of the functions g

d

ℓ

d

m

p

and f , ..., f for F

, as

\+ )<sup>1</sup>

ε/ m /p

1

k

(

ℓ

d]tw (σ, τ) = (1 + ε) · sup{r ∈ R | g(sign(f (π (τ), r, π (σ))), . . .) = 1}.

p

ε

1

ℓ

m

Overall, d]tw corresponds to the ﬁrst value r in R for which the function g of Deﬁnition [10](#br8)

p

ε

outputs a 0. Note that the algorithm also outputs 0 for all larger values in R . Notably, the

ε

function g of Deﬁnition [10](#br8)[ ](#br8)outputs 1 for r/(1 + ε). In the following lemma, we formally show

that d]tw<sub>p</sub>(σ, τ) approximates p-DTW between σ and τ within a factor of 1 + ε.

▶ Lemma 11. Let 0 < ε ≤ 1. For any σ ∈ X<sup>d</sup> and τ ∈ X<sup>d</sup> it holds that

m

ℓ

dtw (σ, τ) < d]tw (σ, τ) ≤ (1 + ε) dtw (σ, τ).

p

p

p

Proof. Let r = d]tw (σ, τ) ∈ R . By deﬁnition the function g of Deﬁnition [10](#br8)[ ](#br8)outputs 1 with

p

ε

σ, τ and r/(1 + ε). Thus dtw<sub>p</sub>(σ, τ) ≤ (1 + ε)r/(1 + ε) = r. As the algorithm outputs 0 for

σ, τ and r it follows that dtw<sub>p</sub>(σ, τ) > r/(1 + ε) implying the claim.

◀

Moreover, from the deﬁnition of d]tw<sub>p</sub>, we conclude that g serves as a membership predicate

for balls deﬁned by d]tw<sub>p</sub>.

▶ Lemma 12. Let ε ∈ (0, 1], τ ∈ X<sup>d</sup> and r ∈ R . For any σ ∈ X<sup>d</sup> the output of the function

ℓ

ε

m

g of Deﬁnition [10](#br8)[ ](#br8)with σ, τ and r corresponds to the decision whether the curve σ is in the

r-ball {x ∈ X<sup>d</sup><sub>m</sub> |

d]tw<sub>p</sub>(x, τ) ≤ r} centered at τ.

Proof. Let r = d]tw (σ, τ) ∈ R . Assume r ≤ r which by Lemma [11](#br8)[ ](#br8)implies that

′

′

p

ε

dtw (σ, τ) ≤ r. Then the function g of Deﬁnition [10](#br8)[ ](#br8)with σ, τ and r outputs 1. Now

p

let r < r ∈ R . However, in this case g with σ, τ and r will by deﬁnition of d]tw

′

ε

p

output 0. Thus membership to a ball range corresponds to the output of the function g of

Deﬁnition [10.](#br8)

◀

We conclude with the main result of this section, namely an upper bound on the VC

dimension of the range space that approximates the p-DTW range space.

▶ Theorem 13. Let ε ∈ (0, 1] and Re<sup>p</sup> = {{x ∈ X<sup>d</sup> | d]tw (x, τ) ≤ r} ⊂ X<sup>d</sup> | τ ∈ X<sup>d</sup>, r > 0}

p

m,ℓ

m

m

ℓ

be the range set consisting of all balls centered at elements of X under d]tw in X . The

d

d

m

p

ℓ

VC dimension of (X , Re<sup>p</sup> ) is at most

d

m

m,ℓ

2(d + 1)ℓ log (12ℓm⌊(m + ℓ)<sup>1</sup><sub>/p</sub>ε<sub>−</sub><sup>1</sup> + 1⌋ + 12m + 12ℓ) = O(dℓ log(ℓmε<sub>−</sub><sup>1</sup>)).

2

Proof. This follows from Theorem [7,](#br6)[ ](#br6)Lemma [9](#br7)[ ](#br7)and Lemma [12,](#br8)[ ](#br8)and the fact that any ball of

radius r > 0 under d]tw coincides with some ball with radius re ∈ R under d]tw . Finally,

p

ε

p

the statement is implied by the injectivity of the functions π and π .

◀

m

ℓ



<a name="br9"></a> 

Conradi, Kolbe, Psarros and Rohde

9

t

y

1

Figure 3 Violated triangle inequality as dtw(s, t) ≈ 12, but dtw(s, x) ≈ 0 (matching in blue),

dtw(y, t) ≈ 0 (red matching) and dtw(x, y) ≈ 3 (green matching).

In this section, we deﬁned a distance function d]tw between curves in X and those in

d

m

p

X<sup>d</sup> that (1 + ε)-approximates dtw and an upper bound on the VC dimension of the range

p

ℓ

space induced by balls of d]tw , thereby producing an approximation of the p-DTW range

space that we make use of below. The bound on the VC dimension is comparable to the one

p

we get for exact p-DTW balls in certain cases (Lemma [43](#br26)[ ](#br26)in Appendix [A).](#br26)[ ](#br26)We emphasize

that the sole purpose of d]tw<sub>p</sub> is to obtain bounds on the size of a sample constituting a

coreset through the knowledge of the VC dimension. At no point do we compute d]tw<sub>p</sub>.

4

Sensitivity bounds and coresets for DTW

To make use of the sensitivity sampling framework for coresets by Feldman and Langberg [[25](#br25)],

we recast the input set T ⊂ X as a set of functions. Consider for any y ∈ X the real

d

m

d

m

valued function f deﬁned on (ﬁnite) subsets of X by f (C) = min dtw (y, c) for C ⊂ X ,

d

d

y

ℓ

y

c∈C

p

ℓ

transforming T into F = {f | τ ∈ T}. To construct a coreset, one draws elements from T

T

τ

according to a ﬁxed probability distribution over T, and reweighs each drawn element. Both

the weight and sampling probability are expressed in terms of the sensitivity of the drawn

element t, which describes the maximum possible relative contribution of t to the cost of any

query evaluation. In our case, as we restrict a solution to a size of k, it turns out that it

suﬃces to analyze the sensitivity with respect to inputs of size k.

▶ Deﬁnition 14 (sensitivity). Let F be a ﬁnite set of functions from P(X<sup>d</sup>) \ {∅} to R. For

ℓ

any f ∈ F deﬁne the sensitivity

f(C)

s(f, F) =

sup

.

~~P~~

P

( )

g C

C={c1,...,c<sub>k</sub>}⊂Z:

g(C)>0 g∈F

g∈F

P

The total sensitivity S(F) of F is deﬁned as

s(f, F).

f∈F

A crucial step in our approach is to show that any (α, β)-approximation for (k, ℓ)-median

under dtw can be used to obtain a bound on the total sensitivity associated to approximate

p

distances. This is facilitated by the following lemma, that is a weaker version of the triangle

inequality, as in general dtw<sub>p</sub> is not a metric (see Figure [3).](#br9)

▶ Lemma 15 (weak triangle inequality [[34](#br26)]). For two curves x and z of complexity m > 0 and

any curve y of complexity ℓ > 0 it holds that dtw (x, z) ≤ m<sup>1</sup> (dtw (x, y) + dtw (y, z)).

/p

p

p

p

![](data:image/jpeg;base64,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)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCACGAAcDASIAAhEBAxEB/8QAFwABAQEBAAAAAAAAAAAAAAAAAAkIA//EACYQAAAFBAIDAAIDAAAAAAAAAAMFBgcIAAECBDl5CRG3ChcSFBX/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8AqKc+ZMBKfkLG0MVE7O6oY5pmOy2QREkG50LGoIcntEIkX6zzcvZyMQAhDlJJVvFqVaGVtXAUj2zkQrvbZxN8tjWVq5fQZiJoecqPy+L4+trorVwofTOeJcKfVIAAjhUueC6MakaEtzjatlfLaPw0u4S2I8N3K1s7F6mNgPXrZy9KDVzkczkTOvSaX3eG9KORzORM69Jpfd4b0oDkczkTOvSaX3eG9KORzORM69Jpfd4b0oDkczkTOvSaX3eG9KORzORM69Jpfd4b0oDkczkTOvSaX3eG9KORzORM69Jpfd4b0oDkczkTOvSaX3eG9KORzORM69Jpfd4b0oDkczkTOvSaX3eG9KORzORM69Jpfd4b0oODiD2F8zcTRLY/xAv49Jo44DXyta1xf3tDi+YImF/WYIuFrXvcMW2OV7Wyva3rG9KhVIvzSmCJ87SQZwihQ6incBlmclBHIlJzpbppvdN1MVeoWydssdpMnKtL9Ut0kZtpBh1UMW/2htjI6FMNH/P3BsgssRFBZRwPEGzq88qZD5MF1tgLXcxYk7ab9dn9jAMFOLUDWKiNJL5I7ZMMV2Cys3uy5SYVGieim4O2MqS/bJwC2+gJlkpSg//Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAbAAcDASIAAhEBAxEB/8QAFwABAQEBAAAAAAAAAAAAAAAAAAkICv/EACAQAAICAwACAwEAAAAAAAAAAAQFBgcCAwgAAQkKFRT/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8Ap83+ZrRF/sKt+MntptZLznFue5zXiWIV0qyP06un12lLP5pvsreWxB1ks4nF68mihaSOLl+QY8zWa/ZmptsLEebCn/G3LA3ze88TlXzzTgM0nXHHZ9tzOU6K/jYjyVWJ6tbmOJDzN82GX5ktpSNGp5M0eD0zdtYflyh2Dju/nZl4bHgRhubqr5Ga8+y9GqSJn6mwoJHKmuzCAqYHUihiRBOc7liWVoL45L2I0TxJxYrp/U1WiFycsjfj+h7HQrG+8WQmjFvOwJpQFNauho50lor9GNeWivZZVuVjieixHpsAkLSMSFlG22IxWpe6FxbQ6PErCG4Rx6LWFuDRlrgmbUY14H//2Q==)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAbAAcDASIAAhEBAxEB/8QAFwABAQEBAAAAAAAAAAAAAAAAAAkICv/EACAQAAICAwACAwEAAAAAAAAAAAQFBgcCAwgAAQkKFRT/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8Ap83+ZrRF/sKt+MntptZLznFue5zXiWIV0qyP06un12lLP5pvsreWxB1ks4nF68mihaSOLl+QY8zWa/ZmptsLEebCn/G3LA3ze88TlXzzTgM0nXHHZ9tzOU6K/jYjyVWJ6tbmOJDzN82GX5ktpSNGp5M0eD0zdtYflyh2Dju/nZl4bHgRhubqr5Ga8+y9GqSJn6mwoJHKmuzCAqYHUihiRBOc7liWVoL45L2I0TxJxYrp/U1WiFycsjfj+h7HQrG+8WQmjFvOwJpQFNauho50lor9GNeWivZZVuVjieixHpsAkLSMSFlG22IxWpe6FxbQ6PErCG4Rx6LWFuDRlrgmbUY14H//2Q==)![ref3]![ref5]![ref6]![ref2]

<a name="br10"></a> 

10

Fast Approximations and Coresets for (k, ℓ)-Median under Dynamic Time Warping

Note that the distance function in question is not dtw<sub>p</sub>, but the (1 + ε)-approximation

d]tw of DTW from before. For any y ∈ X and ε > 0, let fe : P(X ) \ {∅} → R with

d

m

d

ℓ

p

y

e ( ) = min d]tw ( ). Similarly, let e be the set

e

for any

X<sup>d</sup> .

m

f C

y

p y, c

F

{f | τ ∈ T}

T ⊂

c∈C

T

τ

▶ Lemma 16. Let 0 < ε ≤ 1 and let T ⊂ X<sup>d</sup> be the input of size n for (k, ℓ)-median and

m

ˆ

let C = {cˆ , . . . , cˆ } ⊂ X be an (α, β)-approximation to the (k, ℓ)-median problem on T

d

1

ˆ

ˆ

P <sup>k</sup>

ℓ

ˆ

ˆ

ˆ

i

with cost ∆ =

min dtw (τ, cˆ), of size k ≤ βk. For any i ∈ [k] let V = {τ ∈ T |

ˆ

p

τ∈T

cˆ∈C

cˆ∈C

dtw (τ, cˆ ) = min dtw (τ, cˆ)} be the Voronoi region of cˆ , the set of which (breaking ties

p

i

ˆ

p

i

ˆ

P

ˆ

ˆ

i

arbitrarily) partitions T. Let ∆ =

dtw (τ, cˆ ) be the cost of V . For all τ ∈ V let

∈ ˆ

V

p

i

i

i

τ

i



!

ˆ

2α dtw (τ, cˆ )

4

ˆ

8α∆

i

p

i

γ(f<sub>τ</sub>

e ) := ( <sub>)</sub>1/p

mℓ

\+

\+

.

ˆ

ˆ ˆ

∆

∆

|V | |V |

i

i

P

ˆ

Then s(fe , Fe ) ≤ γ(fe ) for any τ ∈ T, and S(Fe ) ≤

γ(fe ) ≤ (mℓ)<sup>1</sup> (4k + 10α).

/p

τ

T

τ

T

τ∈T

τ

ˆ

ˆ

ˆ

i

Proof. Fix an arbitrary set C = {c , . . . , c } ⊆ X , i ∈ [k] and τ ∈ V . Let further

d

ℓ

1

k

i

ˆ

ˆ

i

ε

ˆ

ˆ

B = {σ ∈ V | dtw<sub>p</sub>(σ, cˆ ) ≤ 2∆ /|V |}. Observe that for any τ ∈ B it holds that

i

i

i

i

e ( ˆ) (1 + )2∆

ˆ

ˆ

ˆ

ˆ

f

C ≤

/|V | ≤ /|V |

4∆ . Breaking ties arbitrarily, let ( ) be the nearest

c x ∈ C

τ

i

i

i

i

neighbour of x ∈ X among C.

ˆ

ˆ

P

ˆ

We observe that |B | ≥ |V |/2, since otherwise

dtw (τ, cˆ ) > ∆ . Additionally

i

i

p

i

i

σ

∈ ˆ \ ˆ

V

B

i

P

ˆ

i

note that

fe (C) ≥ ∆/α.

e e

σ

f

∈F

σ

T

ˆ

e

dtw (ˆ

(

))

ˆ

For any σ ∈ B Lemma [15](#br9)[ ](#br9)implies that f (C) ≥ dtw (σ, c(σ)) ≥ p c ,c σ − <sup>4∆</sup>  . Now

i

/p

i

i

σ

p

1

~~ˆ~~

ℓ

|V |

i

as dtw (cˆ , c(σ)) ≥ dtw (cˆ , c(cˆ )) it holds that

p

i

p

i

i





(



!

)

∆ˆ 

~~|~~V ~~|~~

ˆ

2

dtw (ˆ (ˆ )) 4∆

ˆ

∆

ˆ

 X

X

p c , c c

i

i

i

e ( ) max

e ( )

max

i

f

C ≥

f

C ,

≥

−

,

.

ˆ

σ

σ

<sub>ℓ</sub>1<sub>/p</sub>



α



α

|V |

ˆ

i

e e

f

∈F

σ∈B<sub>i</sub>

σ

T

ꢂ

ꢃ

ˆ

2

ˆ

ˆ

ˆ

ˆ

dtw (cˆ ,c(cˆ )) 4∆

∆

, i.e. dtw (cˆ , c(cˆ )) ≥ (2ℓ<sup>1</sup><sub>/p</sub>)

∆+2

α

∆

i

Assume that ~~|~~V ~~|~~

eliminate the dependence on eleme

−

≥

=: δ<sub>i</sub>. To

p

i

/p

i

i

i

p

1

~~ˆ~~

i

i

~~ˆ~~

i

ℓ

|V |

α

α|V |

i

nts in C, we consider the function

2m<sup>1</sup><sub>/p</sub>(ℓ<sup>1</sup><sub>/p</sub> dtw (τ, cˆ ) + x)

p

i

<sup>h</sup>i,τ <sup>: [δ</sup>i<sup>, ∞) → R</sup>><sup>0, x →</sup>

ˆ

ˆ

~~|~~V ~~|~~x 2∆<sub>i</sub>

/p

i

−

2ℓ<sup>1</sup>

and observe that it is monotone and thus its maximum is either h<sub>i,τ</sub> (δ ) or lim

h

(x).

i

x→∞ i,τ

Importantly, as Lemma [15](#br9)[ ](#br9)implies that fe (C) ≤ (1 + ε)m<sup>1</sup> (dtw (τ, cˆ ) + dtw (cˆ , c(cˆ ))),

/p

τ

p

i

p

i

i

~~e~~ <sup>(</sup>

)

we have

f

C

≤ h (dtw (cˆ , c(cˆ ))). If instead dtw (cˆ , c(cˆ )) < δ , then

~~P~~

τ

i,τ

p

i

i

p

i

i

i

e <sup>(</sup>

C

)

f

σ

f

e

σ∈F

e

T

~~e~~ ( )

e ( ) (1 + )

1

(dtw ( ˆ ) +

)

f

C

f

∆

C

ε m <sup>/p</sup> p τ, c <sup>δ</sup>i

τ

τ

i

( )

i

~~P~~

≤ ˆ

<

ˆ

≤ h<sub>i,τ</sub> δ .

e ( )

∆

f

C

/α

/α

e e

σ

f

∈F

σ

T

Thus it folows that

~~e~~ ( )

n

o

f

C

s(f , F

e

e ) =

sup

τ

≤

max ( ) lim ( )

h<sub>i,τ</sub> δ , h<sub>i,τ</sub>

x

~~P~~

τ

T

( )

i

e

x→∞

C={c1,...,c<sub>k</sub>}⊆Z e e f<sub>σ</sub> C

f

σ

∈

F

(



T

!

)

ˆ

2α dtw (τ, cˆ )

4

ˆ

8α∆ 4(mℓ)<sup>1</sup>

/p

= max (mℓ)<sup>1</sup><sub>/p</sub>

p

i

\+

\+

i

,

= γ(fe ).

ˆ

ˆ ˆ

ˆ

i

∆

∆|V |

τ

|V |

|V |

i

i

![ref7]![ref8]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEACEDASIAAhEBAxEB/8QAGAABAQADAAAAAAAAAAAAAAAAAAYECQr/xAApEAAAAwYGAQUBAAAAAAAAAAABAgYAAwQFCRIHCBETFRYUGBkxWZjW/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AOt3ETIKOJKyVC3HOrUAQYqyYx09OksN8x3V0SnjvQAwyxLSHpsdxEnJdo7gfLiLQKUN0dGjZdTIB7AQTwag1TsgnhXBrXebASkLc7KNpC9CG0ofBQ1HQNA1YxgzfbEL9hFT79Yj/AtZYeU9S4erpJrn1vVDFv1OfS2fdQxDzLCpUMpeOiSRPDKxP9LguZkMfZsTKXeXDeVDHO63nd1wGMGw5jGMH//Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEAEIDASIAAhEBAxEB/8QAGAABAAMBAAAAAAAAAAAAAAAAAAUGCQr/xAAsEAAAAgYKAgIDAAAAAAAAAAABAgADBAUGCQgREhMUFhlZmNYHFSFBMlFh/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AOt/yFQEHyTGkTRuNNaYDAZoteTc+zwn42pHZWgmHjrRKYXbCzhya3eodCu3ZVsWLaLJSlC9GqsaQ6pZILnaxLRmCzOyCsZ1ZrCuleJFZKw/EhchDZKH0FY1ftCECQ0xC7hEz7liPQUaYhdwiZ9yxHoKEIDTELuETPuWI9BRpiF3CJn3LEegoQgNMQu4RM+5Yj0FGmIQfgZg8z0Q+wGliIgP8EMg/IIQgaLOyFvWO13u0Iiilv8AXsLIw455PfFPFtwihWoxbwasOrxLa03d81r7sl8vOsWWC2rIEIQP/9k=)![ref6]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEAA4DASIAAhEBAxEB/8QAFgABAQEAAAAAAAAAAAAAAAAAAAcK/8QAJBAAAAYBAwQDAAAAAAAAAAAAAQIDBAUGCAcSEwARFRYJFBf/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A1s6iYBDqLbrXdxzX+QOjDbJCUsJ6lp1kh6xSa8o5OZQYmrQPpj7xEGjybEGH23HGmQheY23uNGxmxeUx/sFjkRyQyr1uCbr0RHkj8hNX/wBJjIMqahH5ndeaevQvjXyx1DN1nHItvaAVHYHbcLp0H//Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEAFQDASIAAhEBAxEB/8QAGAABAQADAAAAAAAAAAAAAAAAAAYFCAr/xAArEAAAAgcIAgIDAAAAAAAAAAABAgADBAUGCRIHCBETGVmY1hQWFSEiYXH/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A637Qbgg2kxtEkbjfVmAQGMXPRtfZ4Ts2vHerQTDp1tJxdsKuH01u+HdBKqVbF5bRSUpQzRwxGIdUskFzuY1ozBZnZBOpKahXevEisuOP0QvoQ4F/WI/1CEDIaYhdwiZ9yxHoKNMQu4RM+5Yj0FCEBpiF3CJn3LEego0xC7hEz7liPQUIQGmIXcImfcsR6CjTELuETPuWI9BQhAaYhdwiZ9yxHoKTzTLYMD5F2Gv/AEzFYoB3Ha81beqraaytDMrys0YEwyBzazK6ftYRWar8cBIQN9LDbGgsVs/YoCC1S2a1fwHi9W0Yztqjn3uPm0Xk1mahZXhEfxbr8ljYRP47uUeGTxmYpVVR8KhIQgf/2Q==)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEABkDASIAAhEBAxEB/8QAGAABAAMBAAAAAAAAAAAAAAAAAAQGCQr/xAAoEAAABQIEBgMBAAAAAAAAAAABAgMEBQYSAAgJEQcTFBYiQRVTkuH/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A63+IWQQeJNaVNW451tQGgzVbJPps9J8Nsx3a1E08dUSmGNpaB7NffERCd9qbLq3FpSlDmjtuNJitMgFY5kqOoNqdkFRskaxPNfYmXcgDaQoUF4lD0G47B7wwwE1LTnNS0tT80TPlqRToMJ2MeKQ9S5oRlYKTIxdJvehlY/sdDrI92ZAqD1tzU+e1UVSvJfcGmPn9qv6/mGGA/9k=)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEAA4DASIAAhEBAxEB/8QAFgABAQEAAAAAAAAAAAAAAAAAAAcK/8QAJBAAAAYBAwQDAAAAAAAAAAAAAQIDBAUGCAcSEwARFRYJFBf/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A1s6iYBDqLbrXdxzX+QOjDbJCUsJ6lp1kh6xSa8o5OZQYmrQPpj7xEGjybEGH23HGmQheY23uNGxmxeUx/sFjkRyQyr1uCbr0RHkj8hNX/wBJjIMqahH5ndeaevQvjXyx1DN1nHItvaAVHYHbcLp0H//Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEAEQDASIAAhEBAxEB/8QAGAABAAMBAAAAAAAAAAAAAAAAAAQGCAr/xAAsEAAAAwUHBAEFAAAAAAAAAAAAAQIDBAUGCQcIERIZWdYTFBaYFRckMUGI/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AOt+0K4IdpM6TNO531qgMhqm2JP0bXKdm147xaSZeW1NKjhsrQHw1++IhDPPlZuXdvGVKUl1TwxOkwumQTWHOTU6gtTtBtHZkvIzvXmhmjFJHlQnwI8qS/BFieBfsAAT9MRO4RU+9sT4CGmIncIqfe2J8BAADTETuEVPvbE+AhpiJ3CKn3tifAQAA0xE7hFT72xPgIaYqdwip8X9YnwEAAa7u52MLsQkqMySu12262PCcYvFUTdbtPv1Bndmh8cYSwKEIj3xUIywVzNzNu4uPafbvD4/NeqvuMEgAB//2Q==)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEAJ4DASIAAhEBAxEB/8QAGAABAQADAAAAAAAAAAAAAAAAAAUGCAr/xAArEAAAAwYFBQACAwAAAAAAAAAAAQIDBQkZWZgEBggR1gcSExUWFBchYXH/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A63+oOgQ+pOdsyZ3PWrEAyGebnpjX2vKfTbUd8tknLq2vas3blVw/G4707oR3drPBfl4jtSlJeU9tzwh1QySbO7BtTiCxO0GtilXYz1Xmhmnff+EJ+CPZP9bn/oAAoSxE1CIn12J8BCWImoRE+uxPgIAASxE1CIn12J8BCWImoRE+uxPgIAASxE1CIn12J8BCWImoRE+uxPgIAASxE1CIn12J8BCWImoRE+uxPgIAASxE1CIn12J8BCWImoRE+uxPgIAASxE1CIn12J8BCWImoRE+uxPgIAASxE1CIn12J8BCWImoRE+uxPgIAASxE1CIn12J8BCWImoRE+uxPgIAAiP2Gu0dmGYNWUQGJo3NtjMJh1Ixeqs27Mkt8QyZmtKfhE7NWff3sl7n2NEpVse2x7J6c9IyNP7/AH8/k6mNXHWn3ro9R6br/wBZv2O4XN4nkzxhvNwu75tzeue7ZTE2DfG+Zt5MK2bsvEXk7iAA/9k=)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEAE0DASIAAhEBAxEB/8QAGAABAAMBAAAAAAAAAAAAAAAAAAUGBwr/xAAtEAAAAwQKAgAHAAAAAAAAAAAAAQIDBQYJBAcIERIUGVmY1hUWEyEjMjNBUf/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDrfrCsCHWTGkTRudtaYDAaoteVOfa4Tq2tHerQTDy2ppUbthZw+m07xDoZ48LOhZukYUpSXxTuvOkOqWSTZ20JqcwWZ2g2lHZqwM7V5oZovL7UJ9CPCkv0V53f0AASGmIncImfcsT6CGmIncImfcsT6CAAGmIncImfcsT6CGmIncImfcsT6CAAGmIncImfcsT6CGmIncImfcsT6CAAGmIncImfcsT6CNVqUslrqZpsZsDtPWvK2Cfa3EyS1rqro97bORDnS91JJwNDhp1+O8kb1V5b5Ns5kHZ+PK/UAA//2Q==)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEACUDASIAAhEBAxEB/8QAGAABAAMBAAAAAAAAAAAAAAAAAAQGCAr/xAAqEAAABAUDAQkBAAAAAAAAAAABAgMFAAQGBxIICRETFBUWFxlBWZjWIv/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDrfuFoEG5NaVNW461twGgzVa5Tz2ek7bajvC1E08dUSmFtpZh8Gz3dDQnninJdrmMSlKHVHjkaS1bZAKt0kqO4NudkFSWSNgnqvwTLyQBxIUKC/koewcjwHvCEBP8ATEL8hG599sR/Aw9MQvyEbn32xH8DCEBqDTXYE9g26t6bNfDUNe0rpUUk5pveoK5XmRULQUrLJSwtjK59ys3YWlQxBm1JPoK5Tiqq/UDPEEIQH//Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEAKkDASIAAhEBAxEB/8QAGgABAAEFAAAAAAAAAAAAAAAAAAIEBQYJCv/EAC0QAAADBgQGAQUBAAAAAAAAAAABAgUJGVmY1gMEBhEHCBIUFiETFWFkcZTi/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AOt/iDyCHxJ1tqTW586rwDQZ6uamdba9J8NuY7xbROnV4vSs2bpVg+G576OyEdXTh5Lu8x0pSkvlPbc8IZTsksZnZPFN4K87Qa8FKujD5rzRhp339IT4Eeyftuf7AAFwhiJmEPPqsTsEIYiZhDz6rE7BAACGImYQ8+qxOwQhiJmEPPqsTsEAAIYiZhDz6rE7BCGImYQ8+qxOwQAAhiJmEPPqsTsEIYiZhDz6rE7BAACGImYQ8+qxOwQhiJmEPPqsTsEAAIYiZhDz6rE7BCGImYQ8+qxOwQAAhiJmEPPqsTsEIYiZhDz6rE7BAACGImYQ8+qxOwQhiJmEPPqsTsEAAIYiZhDz6rE7BEVuxiQhayeDvPd0pUot+bEzLdJGZbl4F7L17IAAUOindasDPMrUeNz2vHGjiMltFnTZLS5nu6YTUSyW0lScg12f4QjvWc0O3SlqZb5sLu04mMXXh9fraF2P5me/o/wAAP/Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEAFcDASIAAhEBAxEB/8QAGAABAAMBAAAAAAAAAAAAAAAAAAUGBwr/xAArEAAAAwUIAgIDAQAAAAAAAAAAAQIDBAUGCQcIERIZWZjWExYUFRgyQVH/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A63rQrgh2kzpMs7nfVqASGqbom/Rtcp2bXjvVpJl1bXKs4bKsB9Nfvp4QjNlZuXy3jKlKS8p4YnSYVTJJtDXJqdQWp2g2juzVkZ3rzQzRiX6oT6EeVJfwsTw/0AASGmIncIqfcsT6CGmIncIqfcsT6CAAGmIncIqfcsT6CGmIncIqfcsT6CAAGmIncIqfcsT6CGmIncIqfcsT6CAAGmIncIqfcsT6CGmIncIqfcsT6CAANbsTuPJsVnp3nn8u78VrJu8MiMNKUbbLfvfJHa/Ys0Mzf20A9ThXkiTlkzw96+WXxmilq8a8cCAAD//Z)![ref7]![ref8]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEACEDASIAAhEBAxEB/8QAGAABAQADAAAAAAAAAAAAAAAAAAYECQr/xAApEAAAAwYGAQUBAAAAAAAAAAABAgYAAwQFCRIHCBETFRYUGBkxWZjW/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AOt3ETIKOJKyVC3HOrUAQYqyYx09OksN8x3V0SnjvQAwyxLSHpsdxEnJdo7gfLiLQKUN0dGjZdTIB7AQTwag1TsgnhXBrXebASkLc7KNpC9CG0ofBQ1HQNA1YxgzfbEL9hFT79Yj/AtZYeU9S4erpJrn1vVDFv1OfS2fdQxDzLCpUMpeOiSRPDKxP9LguZkMfZsTKXeXDeVDHO63nd1wGMGw5jGMH//Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEADYDASIAAhEBAxEB/8QAGAABAAMBAAAAAAAAAAAAAAAAAAUGCAr/xAArEAAAAgcHBAMBAAAAAAAAAAABAgADBAUGCAkHERITFRYZFFmY1hcyQVH/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A637QpBBtJjSJo3GdaoDAZoteTc+zwnZtMdtaCYeOtEphdsLOHZrdpDoV48Kti6towlKUM0brxpDqpkgudrEtGoLU7IKxnVmwK5rxIrJeH1IXYQ4Sh+BeN39QhAkOMQvcIqfeWI+go4xC9wip95Yj6ChCA4xC9wip95Yj6CjjEL3CKn3liPoKEIGjJdJTgl2fURvoJlJs7c9xutldekzF2yfJrlcXStYteow4w7ccmlvRpv6Zra81fnMgApyy3YkIQgf/2Q==)

<a name="br11"></a> 

Conradi, Kolbe, Psarros and Rohde

11

Overall it follows that



!

′

2 dtw ( ˆ )

p τ, c<sub>i</sub>

4

ˆ

8 ∆ˆ

ꢂ

ꢃ

X<sup>k</sup> X

α

α

ˆ

k

S(F<sub>T</sub> ≤ mℓ

e )

(

)<sup>1/p</sup>

\+

\+

i

= (

mℓ<sup>)1/p</sup>

2 + 4 + 8

α

α .

◀

ˆ

ˆ ˆ

∆

∆

|V | |V |

i=1

ˆ

i

i

σ∈V

i

▶ Lemma 17. A weighted ε-coreset for (k, ℓ)-median of T under the approximate distance

d]tw is a weighted 3ε-coreset for (k, ℓ)-median under dtw .

p

p

For T = {τ , . . . , τ } ⊆ X , we write fe = fe .

d

m

1

n

i

τ

j

Proof. Denote by cost(T, C) be the cost associated to dtw , C ⊂ X with |C| = k, and

d

p

ℓ

similarly denote the cost of d]tw by cost(T, C). We will show that

g

p

X

X

(1 − 3ε) cost(T, C) ≤ w(fe) min dtw (τ , c) = w(fe)f (C) ≤ (1 + 3ε) cost(T, C).

i

p

i

i

i

c∈C

τ ∈S

τ ∈S

i

i

P

Denote by cost(T, C) = w(fe) min d]tw (τ , c) the cost associated to the weighted

d

i

c

∈C

p

i

i

ε-coreset for

d]tw . Observe the relations

p

(1 − ε)cost(T, C) ≤ cost(S, C) ≤ (1 + ε)cost(T, C),

g

d

g

T, C ≤

cost(T, C) ≤ cost(T, C) ≤ (1 + ε) cost(T, C),

g

X

X

X

w(f f C ≤ w f f C

e) ( )

( e) e( ) = cdost(

)

(1 + )

ε

w f f C ,

i

( e) ( )

i

i

i

i

i

τ ∈S

τ ∈S

τ ∈S

i

i

i

which together imply

1 + ε

1 − ε

1 + ε

1 − ε

cost(S, C)

X

d

cost(T, C) ≤

cost(T, C) ≤

cost(T, C) ≤

≤

w(fe)f (C)

g

1 − ε

i

i

τ ∈S

i

cost(S, C) ≤ (1 + ε)cost(T, C) ≤ (1 + ε)<sup>2</sup> cost(T, C) ≤ (1 + 3ε) cost(T, C).

d

g

≤

◀

▶ Deﬁnition 18 ([[27](#br25), Deﬁnition 2.3]). Let ε, η ∈ (0, 1) and (X, R) be a range space with ﬁnite

non-empty ground set. An (η, ε)-approximation of (X, R) is a set S ⊆ X, such that for all

R ∈ R

(

ꢄ

ꢄ

~~|~~R~~∩~~X~~|~~ if

ꢄ

ꢄ

ε ·

ε · η,

,

|R ∩ X| ≥ η · |X|

else.

~~|~~R ~~∩~~ X~~|~~ ~~|~~R ~~∩~~ S~~|~~

ꢄ

ꢄ

|X|

−

≤

ꢄ

ꢄ

|X|

|S|

We employ the following theorem for obtaining (ε, η)-approximations.

▶ Theorem 19 ([[27](#br25), Theorem 2.11]). Let (X, R) be a range space with ﬁnite non-empty

ground set and VC dimension D. Also, let ε, δ, η ∈ (0, 1). There is an absolute constant

c ∈ R<sub>>0</sub> such that a sample of

ꢅ

ꢅ ꢆ

ꢅ ꢆꢆ

c

1

1

·

D log

\+ log

η · ε<sup>2</sup>

η

δ

elements drawn independently and uniformly at random with replacement from X is a

(η, ε)-approximation for (X, R) with probability at least 1 − δ.

![ref7]![ref8]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEACEDASIAAhEBAxEB/8QAGQABAAIDAAAAAAAAAAAAAAAAAAIEBggK/8QALBAAAgAEBQIDCQAAAAAAAAAAAQIAAwYRBAUICRIHFhMUMRUXGVdZlpjV1v/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDrd6iaBT1JrKqK3OtXcAoM1ZmOOz16S6b6ju16Jp55oDHLKWyHs3HeyMnTlaXgfN4jiFUeKbRhuXbZAm4DBTDuDbnaF8LIbjL1YFUXlLU8UXsI8VHooubCwvCEBd+GIv1CNz78sT/AxJNsZEdHO4LudzArKxlzdV5eW4Ug8Ji9hDkjWswuLgkXHrCEBvX7rR8w+qP3aP18IQgP/9k=)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEACADASIAAhEBAxEB/8QAFwABAAMAAAAAAAAAAAAAAAAAAAYHCv/EACcQAAAEBQQCAgMAAAAAAAAAAAECAwYEBQgREgAHCRYUIRMVImFx/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/ANb+4NAg7kvZyPca1eQBhi7ppGzs7T22qO6syW6dXE4y1qyHpsd9PKCZYpwXlxGJSlD5RtcYRKuMkFpdBqjyC8nZBOiU2CdV4kTLe/ohehDYv6uP9000Ca8ZIIy+JUDkF5OziQgDipVfmQ35lCxi9CC4e72uHsA1YO3NBQ7WvpvvwK0a+9wzNKPJMiM/dGovtrFcB1UF0fFdDf6dL/t4FMFBUTh/Lh8VSkPmONhaaD//2Q==)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEACADASIAAhEBAxEB/8QAFwABAAMAAAAAAAAAAAAAAAAAAAYHCv/EACcQAAAEBQQCAgMAAAAAAAAAAAECAwYEBQgREgAHCRYUIRMVImFx/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/ANb+4NAg7kvZyPca1eQBhi7ppGzs7T22qO6syW6dXE4y1qyHpsd9PKCZYpwXlxGJSlD5RtcYRKuMkFpdBqjyC8nZBOiU2CdV4kTLe/ohehDYv6uP9000Ca8ZIIy+JUDkF5OziQgDipVfmQ35lCxi9CC4e72uHsA1YO3NBQ7WvpvvwK0a+9wzNKPJMiM/dGovtrFcB1UF0fFdDf6dL/t4FMFBUTh/Lh8VSkPmONhaaD//2Q==)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEADwDASIAAhEBAxEB/8QAGQABAAIDAAAAAAAAAAAAAAAAAAQFBggK/8QAKxAAAAQDBgYCAwAAAAAAAAAAAQIDBQAEBgcICRESExQWGVmY1hUjMkFR/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AOt60K4INpNaVLW431cQCgzVc5zz2ek7Nrx3K1E06dXScW2lWHk2e+HaCatKclxcxpKUobo5ZjhLVhkgs2ySo4guJ2QVJdM2hO9eJEyZh+JC8hDpKH6DMcv7CEBYdMQvcIxPvLEfQYdMQvcIxPvLEfQYQgHTEL3CMT7yxH0GHTEL3CMT7yxH0GEICC44aASMsMwF/wDxM5od1BIEp29XxCACusRIFQTGgy5LIa92XPn9axCHyNpyHYqxC7KNkNKuVMheDvPWmhN1JOvQP1rtq4VpUspxLc0yAtMk6/ANuwypGbRn0JDYPtuLi6TO6PF7aaEB/9k=)![ref6]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEAB8DASIAAhEBAxEB/8QAFgABAQEAAAAAAAAAAAAAAAAAAAYK/8QAKxAAAQIEBQMCBwAAAAAAAAAAAQIFAwQGEQAHCAkWEhMUFSIXISMyQVFh/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/ANb+YWgQ5k1pU1bnWtuA0GqrXKee10nltqO4tRNPLilKi20sw8NnvSGiH19MOS8uY6UpSO6bXMQ1bZIjNslFO4LudoMSXhq6Ieq8ohouPtQngR6Uj8C5t+8MMAdNskQZCPEG4LudrKe17Ymq8rQbxoafcngQva9x8xYgH+YvcvdBHwyrNsrMa0dfeYBp8zoRS2ZuovldGuxdWmdayt7Y+HSHnrb0z6p5sV5UHxHGVkpr6nY7a2GA/9k=)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEAAoDASIAAhEBAxEB/8QAFgABAQEAAAAAAAAAAAAAAAAAAAgK/8QAIxAAAQMCBgMBAAAAAAAAAAAABQEEBgIDAAcIEhMUERUWCf/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDWrmFoGrzGmcwntWtj9A4PXKDXu0iGXeo/5mCxv2zrfWIikf8AjX3pwjJKltjmHcc9W2lNPNc8ecXEEpdR4KIAISJGkBi2AhDB90pE8WQY0sskJmyG2z3y7/h7RJ7w2u08u3r/ABW+TajDAf/Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEAAkDASIAAhEBAxEB/8QAFgABAQEAAAAAAAAAAAAAAAAAAAgK/8QAJRAAAAQEBQUAAAAAAAAAAAAAAgMFBgAEBwgBERITFhQVVZTj/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/ANcVULC+dVCcT9LvPv5ZAnCsKDmEz6fXFcbYqMMZhJpaIgN7h0525vSuQi5VM6w7ZLMMDvi1Z4Uton/Ornv/AChCA//Z)![ref2]

<a name="br12"></a> 

12

Fast Approximations and Coresets for (k, ℓ)-Median under Dynamic Time Warping

▶ Theorem 20. For fe ∈ Fe, let λ(fe) = 2⌈log (γ(f

e<sup>))</sup> , with ( ) the sensitivity bound of

γ fe

⌉

2

ˆ

Lemma [16,](#br10)[ ](#br10)associated to an (α, β)-approximation consisting of k ≤ βk curves, for (k, ℓ)-

P

e<sup>)</sup>

median for curves in X under dtw , Λ =

λ(fe), ψ(fe) = <sub>λ</sub><sup>(</sup><sub>f</sub> and δ, ε ∈ (0, 1). A

d

m

p

e e

Λ

f∈F

sample S of

ꢂ

ꢂ

ꢃꢃ

ˆ

ˆ

Θ ε <sup>2</sup>αk(mℓ)<sup>1</sup><sub>/p</sub> (dℓ log(ℓmε<sub>−</sub><sup>1</sup>))k log(k) log(αn) log(αk(mℓ)<sup>1</sup><sub>/p</sub>) + log(1/δ)

−

elements τ ∈ T, drawn independently with replacement with probability ψ(fe) and weighted

i

i

by w(fe) =

Λ

is a weighted ε-coreset for (k, ℓ)-median clustering of T under dtw with

i

e

p

|S|λ(f )

probability at lea

i

st 1 − δ.

P

Let cost(T, C) =

fe (C) for T = {τ , . . . , τ } ⊆ X , and denote F[e](#br25)[ ](#br25)= {fe , . . . , fe },

g

d

m

τ

1

n

1

n

τ∈T

with fe = fe .

i

τ

j

Our proof relies on the reduction to uniform sampling, introduced by [\[25](#br25)] and improved

by [\[9](#br24)], allowing us to apply Theorem [19.](#br11)[ ](#br11)In the following, we adapt and modify the proof

of Theorem 31 in [[26](#br25)] and combine it with results from [[35](#br26)] to handle the involved scaling,

similarly to Theorem 4 in [\[19\].](#br25)

Proof. The proof of the theorem depends on analyzing diﬀerent estimators for cost(T, C)

g

for C ⊆ X arbitrary with |C| = k. Consider ﬁrst the estimator

d

ℓ

X

X

X

Λ

cost(S, C) = w(fe) · min d]tw (τ , c) = w(fe) · fe(C) =

fe(C)

d

i

p

i

i

i

( )

i

e

c∈C

|S|λ f

τ ∈S

τ ∈S

τ ∈S

i

i

i

i

for cost(T, C). We see that cost(S, C) is unbiased by virtue of

g

d

h

i

|S|

Λ

|S|

( )

f C

X X

X X e

E cost(S, C) =

ψ(fe )

fe (C) =

j

= cost(T, C).

d

g

j

(

)

j

e

|S|

|S|λ f

i=1 τ ∈T

j

i

=1

τ

∈

T

j

j

We next reduce the sensitivity sampling to uniform sampling by letting G be a multiset that

is a copy of Fe, where each f~~e~~∈ Fe is contained |Fe|λ(fe) times and is scaled by

1

, so that

e

<sup>(</sup>e<sup>)</sup>

|F |λ f

~~e~~

<sup>(</sup>~~e~~<sup>)</sup>

|G| = |F|

e Λ and ( e) = ~~|~~F ~~|~~<sub>λ f</sub> . We clearly have

ψ f

|G|

X

X  ~~e~~ ( ~~e~~)

X

~~|~~F~~|~~λ f

g(C) =

f C

e( ) = e( ) = cgost(

f C

T, C .

)

e

( e)

|F|λ f

g∈G

e e

e e

f∈F

f∈F

For a sample S , with |S | = |S|, drawn independently and uniformly at random with

′

′

replacement from G, consider the estimator for cost(T, C) deﬁned by

g

X

~~|~~G~~|~~

~~cost~~(S<sub>′</sub>, C) =

g(C),

|S<sup>′</sup>|

g∈S′

where again C ⊆ X with |C| = k. We see that ~~cost~~(S , C) is unbiased by virtue of

d

ℓ

′

|S<sup>′</sup>|

|S<sup>′</sup>|

ꢇ

ꢈ

~~| |~~ X X

G

f C F λ f

~~| e|~~

e

( )

( )

~~e~~

1

|S<sup>′</sup>|

X X

~~|~~

S

<sup>′</sup>~~|~~

E ~~cost~~(S , C) =

\=

fe(C) = cost(T, C).

′

g

|S<sup>′</sup>|

e

(fe)

|G|

|S |

′

|F|λ

<sup>i=1</sup> e e

<sup>i=1</sup> e e

f∈F

f∈F

![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEABcDASIAAhEBAxEB/8QAFwABAAMAAAAAAAAAAAAAAAAAAAYHCv/EACcQAAAEBAUFAQEAAAAAAAAAAAECBQYDBBESAAcICRQTFRYiQRgk/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/ANb2YWgQcyXm53wOtXcAYYu1Rnlw7Sy21HeLMhvHi2mFNayB4bPdoSCX2w5LlzFpSlDqjSowpK2yAip0lFHcG3OyDEloRrIeq+yGWpAG0hQYXqUPgVGgfcMMBY2XO30XLp5ozy/bO4M+xRu40a2Y+pLyhmqYqCVPJVVdDFmSXNGSCeFQT/6oXGVJSSm/fj9MzDDAf//Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEACgDASIAAhEBAxEB/8QAFwABAQEBAAAAAAAAAAAAAAAAAAUGCv/EACwQAAADBQUHBQAAAAAAAAAAAAECBgADBAUSBwgJERYUFRkhWZjWExcYYXH/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A637Qbgg2krZSLcb6uIAgxV00jZ2dJ2bXjtLIlOne0nGWpWQ6NjtzyglVLuC2uIpKUoeqOWY4iVYZIPpdBvRxBcTsgnclNQ7vXiR2XPPkQughyL9Zj+sYwUOGIXqEYn3diPgLOGIXqEYn3diPgLGMFmzS4d7d2kpRa/My/qudIqZxGgk7RLxGpkUpAgTFMeCVkh0fBb4gJhUUJlD7VDbUV27Ct3TzMYwf/9k=)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEADADASIAAhEBAxEB/8QAFwABAQEBAAAAAAAAAAAAAAAAAAYFCv/EACsQAAADBQUIAwAAAAAAAAAAAAECBgADBAUJBwgREhYTFBUZIVmY1jFhcf/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDrftBuCDaStlItxvq1AEGKumkbOzpOza8dpZEp073KcZalZDo2O4PKCZsruC3uIylKUNqOGIxEqpkg+l0G9GoLU7IJ3JTZHd68SOy449CF0EOBfrEf1jGDQ5Yhe4RU+8sR9BZyxC9wip95Yj6CxjA5Yhe4RU+8sR9BajRdwQ1mC5Qy3d32agi908rpRMTpC0u8jqpDKAjgz0xpYpU/o2A4rKYnoWKg97cbUoAG0L8sYwf/2Q==)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEADEDASIAAhEBAxEB/8QAGQABAAIDAAAAAAAAAAAAAAAAAAIFBgkK/8QALRAAAAQEAwcDBQAAAAAAAAAAAQIDBAAFBgkIERIHExQWGVmYFSHWYWNxktL/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A639oOAQdpNbVJW441bgFBjV00ezs9J7NsR3K1E06dXScZbSsh5NfejygmrSmy4txpKUob0csxwiVWyQWlzNUbgtzsgnRKbQnivEiZc8/YheQhyL9Mx/MIQFh0xC9wi595Yj8Bh0xC9wi595Yj8BhCAdMQvcIufeWI/AYi3tqCwm8uXC4BczeEaqFf8I/xV8SzcnZrJKFQdoDQZN81WABScoiYu9ROcmoueYIQGxr1Nz9v9B/qEIQH//Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEACQDASIAAhEBAxEB/8QAGQABAAIDAAAAAAAAAAAAAAAAAAIFBggK/8QAKxAAAQIEBQMCBwAAAAAAAAAAAQIDAAQFBgcICRESFBUWIVETGTJBWZjW/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AOt7ELIIcSb0uW9znV1ALDVd1Tnq2u08Nsx3i1k26t3is021aD4bPdnpCOXFuS6uY4pSkfFO25wmlaZIepsk6dQXU7QXJdtXBvNeUNo3H0oT4EeKR9hudveEICw+WIn8hGp9+2J/gYg5pjhDbixqD6npKEKUAc2JIJSkkbjwL1G49R7QhAbj5fcMHMH8MqXYxxKxXxTVIVCsTirzxkvDzi/agqp1B6dUxUri7dTOqlpEu9LTmujb6aTbaY5L4cyhCA//2Q==)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEACYDASIAAhEBAxEB/8QAGAABAAMBAAAAAAAAAAAAAAAAAAUGBwr/xAArEAAABAQFAgYDAAAAAAAAAAABAgMFAAQGCAcJERIWExQVFxlZmNYyQVH/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A638QrBBxJrSpq3G9bMBoM1WuU89npPDa47i1E08dUSmFtpZh4bPeENCe/anJd3MbSlKHVHTUaQ1ZZILNskqOYLmdkFSXTNsTuvEiZNQ/EheBDtKH6DUdP7CEBIemIX3CMz75Yj9Bh6YhfcIzPvliP0GEIDabe7U1bdXx8qktzV2+OJ6gbFWAzFcTjN5mU20FRckZwjuxNfG2QW99EsoEkZw7hbWQXmZfoh1d5UIQH//Z)![ref6]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEAC4DASIAAhEBAxEB/8QAFwABAQEBAAAAAAAAAAAAAAAAAAYFCv/EACsQAAADBQcEAgMAAAAAAAAAAAECBQADBAYHCAkREhMUFhlZmNYVGCFhcf/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDrfqDYEGpM7TJO421bwCQxm5UjVs8p02tHcWkmXTvcpxTZVQeGx3w6QTNldwW7iMpSlDVHDEYhKuyQfJ0G9G8FvOyCdyU2R3avEjsuOP4IXgQ4F/WI/wBYxg0OmIXuEXn3liPoLOmIXuEXn3liPoLGMDpiF7hF595Yj6C11SWxaejk8qSz9t7bVVCRMsxCORFrNXnnKDBb1RSVA6rApwyombdZcimBBw0drH0oGOUYfSNuc7sxg//Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEAC4DASIAAhEBAxEB/8QAFwABAQEBAAAAAAAAAAAAAAAAAAYFCv/EACsQAAADBQcEAgMAAAAAAAAAAAECBQADBAYHCAkREhMUFhlZmNYVGCFhcf/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDrfqDYEGpM7TJO421bwCQxm5UjVs8p02tHcWkmXTvcpxTZVQeGx3w6QTNldwW7iMpSlDVHDEYhKuyQfJ0G9G8FvOyCdyU2R3avEjsuOP4IXgQ4F/WI/wBYxg0OmIXuEXn3liPoLOmIXuEXn3liPoLGMDpiF7hF595Yj6C11SWxaejk8qSz9t7bVVCRMsxCORFrNXnnKDBb1RSVA6rApwyombdZcimBBw0drH0oGOUYfSNuc7sxg//Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEABcDASIAAhEBAxEB/8QAFwABAAMAAAAAAAAAAAAAAAAAAAYHCv/EACcQAAECBAYCAgMAAAAAAAAAAAECBgMEBQcACAkSExYRFBUYISJB/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/ANbNw8gZuK73U9jnX1AmKXZUam4FtO3WY/rDJb0SZUYppTXoPTJ74miQeTjgSBmpjjhIho5VbfOIrJaY4iScpEOoNqeJK5WXWUozYFKE7oSDtSnoX6pT58JH8AAwwwFzWRyQJsm+YT4+3WeC7ZhUqo0sNK9t/e+shfyAgg1FdALTpW6qSXCDITftj1zEinjXv/DDDAf/2Q==)

<a name="br13"></a> 

Conradi, Kolbe, Psarros and Rohde

13

ꢉ

ꢄ

ꢊ

ꢄ

We now assume that S<sub>′</sub> =

1

· feꢄτ ∈ S , which yields

i

ꢄ <sup>i</sup>

e

<sup>(</sup>e<sup>)</sup>

|F |λ f

i

e Λ X

X

Λ

|F|

~~cost~~(S<sub>′</sub>, C) =

g(C) =

fe(C) = cost(S, C).

(I)

d

(fe)

i

|S<sup>′</sup>|

|S|λ

g∈S′

τ ∈S

i

i

For H ⊆ G, C ⊆ X with |C| = k and r ∈ R , we deﬁne range(H, C, r) = range(G, C, r)∩

d

ℓ

≥0

H = {g ∈ H | g(C) ≥ r}. For all such C ⊆ Z and all H ⊆ G, we have that

Z

Z

X

X

∞

(g(C) ≥ r) dr = <sub>0</sub><sup>∞</sup>|range(H, C, r)| dr,

g(C) =

(II)

0

g∈H

g∈H

where all of the involved functions are integrable. Consider now the range space (G, R)

over G, where R = {range(G, C, r) | r ∈ R , C ⊆ Z, |C| = k}. For the following, we

≥0

apply Theorem [19](#br11)[ ](#br11)with the given δ, ε/2 and η = 1/Λ, so as to guarantee that S is a

′

(1/Λ, ε/2)-approximation of (G, R). Given C ⊆ Z with |C| = k, we compute that

ꢄ

ꢄ

ꢄ

ꢄ

ꢄ

ꢄ

ꢄ

ꢄ

ꢄX

X

ꢄ

ꢄ

ꢄ [<sup>(I)</sup>](#br13)[</sup> ](#br13)ꢄ

ꢄ

~~|~~G~~|~~

cost(T, C) − cost(S, C) [=](#br13)[ ](#br13)cost(T, C) − ~~cost~~(S , C) = ꢄ g(C) −

g(C)ꢄ

g

d

g

′

ꢄ

ꢄ

ꢄ

ꢄ

ꢄ

ꢄ

|S<sup>′</sup>|

~~ꢄ~~

ꢄ

g∈G

g∈S′

ꢄZ

Z

ꢄ

ꢄ

∞

∞

ꢄ

~~|~~G~~|~~

= ꢄ |range(G, C, r)| dr −

|range(S , C, r)| drꢄ

′

ꢄ

ꢄ

|S<sup>′</sup>|

0

0

ꢄZ

ꢄ

ꢄ

∞

ꢄ

~~|~~G~~|~~

= ꢄ |range(G, C, r)| − |range(S , C, r)| drꢄ

′

ꢄ

ꢄ

|S<sup>′</sup>|

~~|~~G~~|~~

|S<sup>′</sup>|

0

Z

ꢄ

ꢄ

∞ ꢄ

ꢄ

ꢄ range(

)

range(

) ꢄ d

′

≤

|

G, C, r | −

|

S , C, r | r.

ꢄ

ꢄ

0

The monotonicity of |range(G, C, r)| implies that R (C) = {r ∈ R | |range(G, C, r)| ≥

1

≥0

η · |G|} and R<sub>2</sub>(C) = R <sub>0</sub> \ R<sub>1</sub>(C) are intervals. Denoting r (C) = max g(C), we have that

≥

u

g∈G

for r ∈ (r (C), ∞), it holds that |range(G, C, r)| = 0. Therefore,

u

ꢄ

ꢄ

ꢄ

′

ꢄ

~~|~~

range(G, C, r)~~|~~ ~~|~~range(S , C, r)~~|~~ ε ~~|~~range(G, C, r)~~|~~

ꢄ

ꢄ

−

≤

,

ꢄ

ꢄ

2

|G|

|S<sup>′</sup>|

|G|

for r ∈ R (C), since S is a (1/Λ, ε/2)-approximation for (G, R) and similarly for r ∈ R (C).

′

1

2

Thus,

Z

ꢄ

ꢄ

Z

ꢄ

ꢄ

ꢄ

ꢄ

ꢄ

ꢄ

~~|~~G~~|~~

ε

2

cost(T, C) − cost(S, C) ≤ ꢄ|range(G, C, r)| − |range(S , C, r)|ꢄ dr +

η|G|dr

g

d

′

ꢄ

ꢄ

ꢄ

ꢄ

|S<sup>′</sup>|

R1

R2

(C)

Z∞

<sup>r</sup>Z

u

ε

2

εη|G|

≤

|range(G, C, r)| dr +

dr

2

0

0

X

C

( )

ε

2

<sup>εη|G|r</sup>u

\=

g(C) +

,

(III)

2

g∈G

where the last equality is due to [(II)](#br13). We now bound the last term in [(III)](#br13)[ ](#br13)with the help

of the sensitivity bounds derived in Lemma [16,](#br10)[ ](#br10)where we use that γ(fe) ≤ λ(fe). For each

g ∈ G, we have

1

~~e~~( )

g(C)

f C

1

1

~~e~~

<sup>(</sup>~~e~~<sup>)</sup>

= |F |λ f

≤

λ(fe) =

,

~~P~~

~~P~~

g(C)

e(C) |Fe|λ(fe)

|Fe|

f

g∈G

e e

f∈F

![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEACkDASIAAhEBAxEB/8QAGQABAAIDAAAAAAAAAAAAAAAAAAEFBgkK/8QALxAAAAMGBAMHBQAAAAAAAAAAAQIEAAMFBgkRCBIUFgcTIRUZJjFZYWRxlJjW4v/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDrf4g4BB4kztMk7jjVqASGM3RRbGzynw2xHbWkmXTvcpxhsqwHZq7seEEzZXaLVqMpSlDmja44RCqZIPocjejUFqdkE7kpsjvFeJHZb36ELsIbF9rj9WMYLDuxC+oRU+/LEf0FoNTFKACIVB6n1wARC+LEbdAv18AsYwVEv04RMqdxJ9j1qTKTw+KPxBKoxR8xErIhVOT8hcn2MGocLMpSriZyagoWuTzbaVofmLvuP4Yxg//Z)![ref6]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEAB4DASIAAhEBAxEB/8QAGAABAAMBAAAAAAAAAAAAAAAAAAQGCAr/xAAlEAAABAUEAwEBAAAAAAAAAAABAgMEBQYIERIABwkWExQVFzH/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A63dxKBR3JnOaJ3GtXkAkMZsiD6OnlLbao7q0kS8dYAMMMleAdNffIg6eVkmXtuMSlKHlG2qZDOMgFYcyVHkF5OyCo2RPgnVgJUy5EAcSF6EOJS/woXGwAAXHTTQQo7xrqQxs3VS5AeTRcVnzNscjuqsV0wIu5RSFQpehlsqnn5Ej3HBQpDWG1h01TlSaSn2LTFFS1JVY719jhjRh82oHeP8ASIVAQZPl3AvJbZ9cg3y3zwygpvnHkX9hECkwJa4tNB//2Q==)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEADADASIAAhEBAxEB/8QAFwABAQEBAAAAAAAAAAAAAAAAAAYFCv/EACsQAAADBQUIAwAAAAAAAAAAAAECBgADBAUJBwgREhYTFBUZIVmY1jFhcf/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDrftBuCDaStlItxvq1AEGKumkbOzpOza8dpZEp073KcZalZDo2O4PKCZsruC3uIylKUNqOGIxEqpkg+l0G9GoLU7IJ3JTZHd68SOy449CF0EOBfrEf1jGDQ5Yhe4RU+8sR9BZyxC9wip95Yj6CxjA5Yhe4RU+8sR9BajRdwQ1mC5Qy3d32agi908rpRMTpC0u8jqpDKAjgz0xpYpU/o2A4rKYnoWKg97cbUoAG0L8sYwf/2Q==)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAMAAgDASIAAhEBAxEB/8QAFgABAQEAAAAAAAAAAAAAAAAAAAcK/8QAIBAAAQQCAwADAAAAAAAAAAAABQIEBgcBAwAICRESFP/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDVNStqWJJ/TvunVJmbyZ7XFWVVQL+Hw1RBaI4FczECzOnzf5E42oUZIu9GGyt6fjLge6eYznXlf145YYD53UdW3bCbdyozMOwWLZsV4+czIOTvmwylWHG28MSAigzysHZVcT2gIgPJ7MQkMtipjFN7Ue4D6m+xi3yhwP/Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEABcDASIAAhEBAxEB/8QAFwABAAMAAAAAAAAAAAAAAAAAAAYHCv/EACcQAAECBAYCAgMAAAAAAAAAAAECBgMEBQcACAkSExYRFBUYISJB/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/ANbNw8gZuK73U9jnX1AmKXZUam4FtO3WY/rDJb0SZUYppTXoPTJ74miQeTjgSBmpjjhIho5VbfOIrJaY4iScpEOoNqeJK5WXWUozYFKE7oSDtSnoX6pT58JH8AAwwwFzWRyQJsm+YT4+3WeC7ZhUqo0sNK9t/e+shfyAgg1FdALTpW6qSXCDITftj1zEinjXv/DDDAf/2Q==)![ref9]![ref9]![ref9]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEACIDASIAAhEBAxEB/8QAGAABAAMBAAAAAAAAAAAAAAAAAAQGCQr/xAAsEAAAAwYDBwUBAAAAAAAAAAABAgUAAwQGEhQICREHExUWGSFBIjFRWZjW/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AOt7aFgFHaTOczzuONXMAkMZtUY5cPKezbEdytJMvHfCUwpsrIPJsdwhId1UuoK7iKClKG9HTVqUl5ZAPU6CeDmDZnZBeQzk9DvFhS7LUQo0kKEhekoexQ1HQO2rGME/piF+wjM+/WI/wLOmIXzmD5noh5AcWIiAh8CHIPcB8sYwaJI0pcGR0pICZpuVASk2BTQU1lZvldRsYV1C36rG2zq8UovdXEdFbp3cRTx69oJXSBjGD//Z)![ref9]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEAEMDASIAAhEBAxEB/8QAGQABAAIDAAAAAAAAAAAAAAAAAAMFBggK/8QALBAAAQEGAwgCAwEAAAAAAAAAAgEAAwQFBgkHCBESExQWGVmY1hUhIkFhcf/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDrfxByCLiTW1SVuudW4BQa1dNI2dnSeG2Y7laiadN7smstpWQ8mx3w8oDa2XcFxcRsiIpvV01XCJVbJR9LoN6twW52Cm5Eth3mvUHY66/QDyEug/zVf9YxgsOmIPcIufeWK+gs6Yg9wi595Yr6CxjA6Yg9wi595Yr6CzpiD3CLn3livoLGMDpiD3CLn3livoLRvrZCOnT16Nwe56pO3ZmKFmxVRVQFSRCTkL7TVPtP2jGMGyGFGAx4e0DIaQLGzMBXRSn5UiqvEbEdKnrObFMZ1MZqRTmefDQXHLBrHLL4BeFdcPK4SChPz3G8MxjB/9k=)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEAEQDASIAAhEBAxEB/8QAGAABAAMBAAAAAAAAAAAAAAAAAAQGCAr/xAAsEAAAAwUHBAEFAAAAAAAAAAAAAQIDBAUGCQcIERIZWdYTFBaYFRckMUGI/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AOt+0K4IdpM6TNO531qgMhqm2JP0bXKdm147xaSZeW1NKjhsrQHw1++IhDPPlZuXdvGVKUl1TwxOkwumQTWHOTU6gtTtBtHZkvIzvXmhmjFJHlQnwI8qS/BFieBfsAAT9MRO4RU+9sT4CGmIncIqfe2J8BAADTETuEVPvbE+AhpiJ3CKn3tifAQAA0xE7hFT72xPgIaYqdwip8X9YnwEAAa7u52MLsQkqMySu12262PCcYvFUTdbtPv1Bndmh8cYSwKEIj3xUIywVzNzNu4uPafbvD4/NeqvuMEgAB//2Q==)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEACYDASIAAhEBAxEB/8QAGAABAAMBAAAAAAAAAAAAAAAAAAUGBwr/xAArEAAABAQFAgYDAAAAAAAAAAABAgMFAAQGCAcJERIWExQVFxlZmNYyQVH/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A638QrBBxJrSpq3G9bMBoM1WuU89npPDa47i1E08dUSmFtpZh4bPeENCe/anJd3MbSlKHVHTUaQ1ZZILNskqOYLmdkFSXTNsTuvEiZNQ/EheBDtKH6DUdP7CEBIemIX3CMz75Yj9Bh6YhfcIzPvliP0GEIDabe7U1bdXx8qktzV2+OJ6gbFWAzFcTjN5mU20FRckZwjuxNfG2QW99EsoEkZw7hbWQXmZfoh1d5UIQH//Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEAC4DASIAAhEBAxEB/8QAFwABAQEBAAAAAAAAAAAAAAAAAAYFCv/EACsQAAADBQcEAgMAAAAAAAAAAAECBQADBAYHCAkREhMUFhlZmNYVGCFhcf/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDrfqDYEGpM7TJO421bwCQxm5UjVs8p02tHcWkmXTvcpxTZVQeGx3w6QTNldwW7iMpSlDVHDEYhKuyQfJ0G9G8FvOyCdyU2R3avEjsuOP4IXgQ4F/WI/wBYxg0OmIXuEXn3liPoLOmIXuEXn3liPoLGMDpiF7hF595Yj6C11SWxaejk8qSz9t7bVVCRMsxCORFrNXnnKDBb1RSVA6rApwyombdZcimBBw0drH0oGOUYfSNuc7sxg//Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEABQDASIAAhEBAxEB/8QAFwABAAMAAAAAAAAAAAAAAAAAAAcICv/EACUQAAAEBQQDAQEAAAAAAAAAAAECBAUAAwYHEQgTFiEJEhQVGf/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDW7cTQKNyayqitx1q+QCgxqxxXPp6StvqO4vRNPHmgBhbKWYeGrvyGcntiWh+tR6gUobo4iPEPjOETITf0D8m4ZZQU4DVbgoCU8jEkA4F0nHPcvI5wHfUIQFjbX6VzWkp9VTpNSOq65H3u09+M+3YvDzGo0xlaJvRflI3TjrdssqYrcVQkQbB9lWsXzt430epEIQH/2Q==)![ref2]

<a name="br14"></a> 

14

Fast Approximations and Coresets for (k, ℓ)-Median under Dynamic Time Warping

P

where fe∈ Fe is the function that g is a copy of, implying that r (C) ≤

1

g(C). Thus,

u

e

g∈G

|F |

εη|G|r (C) ε 1

1

X

ε

2

X

u

e Λ

( ) =

g C

( )

g C .

≤

|F|

2

2 Λ

e

|F|

g∈G

g∈G

P

All in all, [(III)](#br13)[ ](#br13)implies that |cgost(T, C) − cdost(S, C)| ≤ ε ·

g(C) = ε · cgost(T, C) for all

g∈G

C ⊆ Z with |C| = k, with probability at least 1 − δ, so S is an ε-coreset for the approximate

distance function.

By Lemma [17,](#br11)[ ](#br11)upon rescaling ε by 1/3, it remains to show the asserted bounds on the size of

. By the use of Theorem [19,](#br11)[ ](#br11)we need a bound on both Λ and the VC dimension of the range

′

S

ˆ

space (G, R). We observe that γ(fe) ≤ λ(fe) ≤ 2 · γ(fe), so that Λ ≤ 2Γ(Fe) = O(αk(mℓ)<sup>1</sup><sub>/p</sub>),

by Lemma [16.](#br10)[ ](#br10)Moreover, the VC dimension of (G, R) lies in O(dℓ log(ℓmε <sup>1</sup>)k log(k) log(αn))

−

by Lemma [21](#br14)[ ](#br14)below, concluding the proof.

◀

▶ Lemma 21. The VC dimension of the range space (G, R) lies in O(Dk log(k) log(αn)),

with D = O(dℓ log(ℓmε <sup>1</sup>)) the VC dimension of (X , Re<sup>p</sup> ).

−

d

m

mℓ

The proof is an adaptation of the methods of [\[35,](#br26)[ ](#br26)Lemma 11] and [\[19,](#br25)[ ](#br25)Lemma 2].

Proof. We ﬁrst consider the simple case that for all fe∈ Fe, λ(fe) = ec, so that the scaling of

the elements of G is uniform and can be ignored in the context of the VC dimension. For

r ≥ 0, c ∈ X<sup>d</sup>, let B

e

( ) =

c

{σ ∈

X<sup>d</sup> d]tw (

|

)

. The range space ( ) can

G, R

p σ, c ≤ r} ∩ X

ℓ

r,X

m

S

then alternatively be described as (T, {T \

Be (c))|C ⊂ X , |C| = k, r ∈ R }), which

d

≥0

c∈C r,T

<sup>ℓ</sup>S

in turn has VC dimension at most equal to that of (X , {X

\

Be |Be , ..., Be }), with

d

m

d

m

i∈

[k]

i

1

k

each Be of the form Be

d

(c) for a c ∈ X . The last range space has the same VC as its

d

ℓ

i

r,X

complementary range spa

m

ce, which by the k-fold union theorem has VC dimension at most

2Dk log (3k) ≤ cDk log k ∈ O(Dk log k) [\[5,](#br24)[ ](#br24)Lemma 3.2.3].

2

Let now t denote the number of distinct values {c , ..., c } of λ(fe), as fe ranges over Fe

1

t

and partition G into the sets {G , ..., G } such that for all g ∈ G there is a fe ∈ Fe with

1

t

i

g =

1

f

e =

1

f

e. Assume, for the sake of contradiction, that

G

<sup>′</sup> ⊂

G

is a set with

e

<sup>(</sup>e<sup>)</sup>

e

|F |λ f

|F |c

i

|G<sup>′</sup>| > t · cDk log k that is shattered by R. Consider the sets G′ = G<sup>′</sup> ∩ G as well as induced

i

i

range spaces R = G ∩R on each G for i ∈ [t]. Since the G are disjoint, each G is shattered

′

i

i

i

i

i

by R and there must exist at least one j ∈ [t] such that |G | ≥ ~~|~~G<sup>′</sup>~~|~~ > <sup>t</sup>·cD<sub>k</sub> <sup>log</sup> <sub>k</sub> = cDk log k.

′

j

i

t

t

However, this contradicts the VC dimension of (G , R ) in the case that the scaling of the

j

j

functions in G is uniform, established above. We now derive explicit bounds on t. By

ˆ

Lemma [16,](#br10)[ ](#br10)for τ ∈ V with i ∈ [k ], we have

′

i

ꢅ

ꢆ

4

4

ˆ

8α

(mℓ)<sup>1</sup><sub>/p</sub> ≤ γ(fe) ≤ (mℓ)<sup>1</sup> 2α +

\+

,

/p

ˆ

ˆ

|V |

|V | |V |

i

i

i

implying that the number of distinct values of λ(fe) is that number for ⌈log (γ(fe))⌉, which is

2

ꢅ

ꢅ

ꢆꢆ

ꢅ

ꢆ

4

ˆ

8α

ˆ

4

ˆ

log (mℓ)<sup>1</sup> 2α +

\+

<sup>− log</sup>2 (mℓ)<sup>1</sup>

/p

/p

2

|V | |V |

|V |

i

i

i



!

ꢅ

ꢆ ꢅ

ꢆ

−

1

ꢂ

ꢃ

8α

4

ˆ

4

ˆ

nα

2

= log<sub>2</sub> 2α +

\+

≤ log<sub>2</sub>

\+ 2α + 1 .

ˆ

|V | |V | |V |

i

i

i

Thus, the VC dimension of (G, R) is at most 2Dk log (3k) log ( + 2α + 1). Theorem [13](#br8)

nα

2

2

2

concludes the proof.

◀

We remark that in the limit p → ∞, the constructed coreset has a very similar size as a

recent construction for coresets for the Fréchet distance [\[19\].](#br25)

![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEABEDASIAAhEBAxEB/8QAFwABAAMAAAAAAAAAAAAAAAAAAAYJCv/EACcQAAAEBAYBBQAAAAAAAAAAAAECAwQABQYSBwgJERMWFRQYY5LS/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/ANb2IWQQcSa0qWtxzq6gFBmq6Zvp2ek8Nsx3VqJp06tpxltKyHpr7w8oJdamy9W4tKUoco7bjB5FpsGaqSF57/dTByVFNKYixdZqeaXrHaKoqlbLtuiF5GawXJOELw5ETnJeXfeEICyrybn4/oP6hCEB/9k=)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEAEMDASIAAhEBAxEB/8QAGQABAAIDAAAAAAAAAAAAAAAAAAMFBggK/8QALBAAAQEGAwgCAwEAAAAAAAAAAgEAAwQFBgkHCBESExQWGVmY1hUhIkFhcf/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDrfxByCLiTW1SVuudW4BQa1dNI2dnSeG2Y7laiadN7smstpWQ8mx3w8oDa2XcFxcRsiIpvV01XCJVbJR9LoN6twW52Cm5Eth3mvUHY66/QDyEug/zVf9YxgsOmIPcIufeWK+gs6Yg9wi595Yr6CxjA6Yg9wi595Yr6CzpiD3CLn3livoLGMDpiD3CLn3livoLRvrZCOnT16Nwe56pO3ZmKFmxVRVQFSRCTkL7TVPtP2jGMGyGFGAx4e0DIaQLGzMBXRSn5UiqvEbEdKnrObFMZ1MZqRTmefDQXHLBrHLL4BeFdcPK4SChPz3G8MxjB/9k=)![ref9]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEAAwDASIAAhEBAxEB/8QAFgABAQEAAAAAAAAAAAAAAAAAAAgK/8QAJhAAAAMHBAIDAAAAAAAAAAAAAgUGAQMEBwgREhMVFhcAFAkhJP/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDW9MKgRsyVopVu2tX5AEGJXGccdjSctqjuLIlOje4jaWpUh4bHbOUAyxdwXtxGIQhZqttdtSU5Se6PlgWoftKcc4LR8cccwnmt+wF0Lc9EW1iUG2FNycv0bFkD6n5WPXzNV5n9PHgf/9k=)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEABQDASIAAhEBAxEB/8QAFwABAAMAAAAAAAAAAAAAAAAAAAcICv/EACcQAAAEBAUFAQEAAAAAAAAAAAECBAUAAwYHCBESFiEJExQVIhkj/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/ANblxcAo3JrKqa3HGr1AKDGrHBe+npK2+I7a9EU+eaAGFspdg2au9SzE1ZS0PlqNIFKHdHTnEeoumdnMQm/QPqb8sZFWkMVvyAgeT/EA2FwnHPmXnzkH1xCEBZi1GGI1nacV00TERijud7B4nv56gu/dretTJjLEDcg9QjddvtnYY0xWwqlI3+ObsrFrhP7xvJ0EQhAf/9k=)![ref9]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEACYDASIAAhEBAxEB/8QAGAABAAMBAAAAAAAAAAAAAAAAAAUGBwr/xAArEAAABAQFAgYDAAAAAAAAAAABAgMFAAQGCAcJERIWExQVFxlZmNYyQVH/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A638QrBBxJrSpq3G9bMBoM1WuU89npPDa47i1E08dUSmFtpZh4bPeENCe/anJd3MbSlKHVHTUaQ1ZZILNskqOYLmdkFSXTNsTuvEiZNQ/EheBDtKH6DUdP7CEBIemIX3CMz75Yj9Bh6YhfcIzPvliP0GEIDabe7U1bdXx8qktzV2+OJ6gbFWAzFcTjN5mU20FRckZwjuxNfG2QW99EsoEkZw7hbWQXmZfoh1d5UIQH//Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEABoDASIAAhEBAxEB/8QAGAABAAMBAAAAAAAAAAAAAAAAAAQGCQr/xAAmEAAABQIFBQEBAAAAAAAAAAABAgMFBgQSAAgJERMHFBUWQSEi/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AOt7qFkEHqTM5POBzq6gEDGWuNc+HiXTbMd6tCI8dW0wtsWYPTa7xDQS+1Oi7uotKUoco7bjSmrTIBVuolR1BtTsgqUyRrE819iZdyANpChAv5KHwNx2D7hhgJNTpk8FNUrF1BtTwxkqdc5SnzYXEESonEAMUYF+hv8ANwxorBo+aPQmHsBn6RyAzHFo+zmfpM6eVkj2ZsaaSiF3kDpwIeSe3IUBrXWv4Ee8r1qio4k+SwrDAf/Z)![ref6]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEADQDASIAAhEBAxEB/8QAGQABAAIDAAAAAAAAAAAAAAAAAAQFBggK/8QALBAAAAIIBAUFAQAAAAAAAAAAAQIAAwQFBgcJEggREyEUGVmY1hUWFyJBI//EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDremFgEGZMZxPHA41agEBjFrxbn4eEpbYjva0EQ8dbaYXbCzg9mt3pDoJfarYuLaLSlKGqOWY4U6qZALXcxLRqDVOyCsZlRrFeK+xWXMgDaQoQF9Sh+BmOQfqEIE/liF6hFT7uxHwFHLEL1CKn3diPgKEIDliF6hFT7uxHwFKs1M3SeTSAVA6nB9F3KF4amK2+8ddo/ktzgIL1P0zs23MbffYhA3CkjJlbIGD2iCFM354ThK0vxsiE0WT3j75DjZUd4MjvZBdCl/ekuexxMYO4q9hd/CG4dpbG9brH4i0hCED/2Q==)![ref8]![ref8]![ref8]![ref8]![ref8]![ref8]![ref8]![ref8]![ref8]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEABMDASIAAhEBAxEB/8QAFwABAAMAAAAAAAAAAAAAAAAAAAcICv/EACYQAAAFAgUFAQEAAAAAAAAAAAECAwQGBQcACBETFgkSFBUhGCP/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A1u3EyCjcmZSibjnV6gEDGWVF9XTxK2+Y7i8Jjx1QAw0yLUHhr71FHJ3aJsfLcdoFKG6OmIwadNk6rhsyPn+6mAo+lbPQH9U/0Kcy7RPZKfgnxsUVdwqWnxVNM/cPZoLDAX+tPZ8LWW+jsCC6F4LihHk6glzO6c25dPq4L6rv6oKshkXrKf7JZoL4aeyP4aOxTGjJpofx9w7DDAf/2Q==)

<a name="br15"></a> 

Conradi, Kolbe, Psarros and Rohde

15

5

Linear time (O((mℓ)<sup>1</sup><sub>/p</sub>), 1)-approximation algorithm for (k, ℓ)-median

In this section, we develop approximation algorithms for (k, ℓ)-median for a set T ⊂ X<sub>d</sub>

m

of n curves. For this, we approximate DTW on T by a metric using a new inequality for

DTW (Lemma [22).](#br15)[ ](#br15)This allows the use of any approximation algorithm for k-median in

metric spaces, leading to a ﬁrst approximation algorithm of the original problem. However,

computing the whole metric space would take O(n<sup>3</sup>) time. We circumvent this by in turn using

the DTW distance to approximate the metric space. Combined with a k-median algorithm

in metric spaces [\[29\],](#br25)[ ](#br25)we obtain a linear time (O((mℓ)<sup>1</sup><sub>/p</sub>), 1)-approximation algorithm.

5\.1 Dynamic time warping approximating metric

We begin with the following more general triangle inequality for dtw , which motivates

p

analysing the metric closure of the input set. While dtw<sub>p</sub> does not satisfy the triangle

inequality (see Figure [3),](#br9)[ ](#br9)the inequality shows it is never ‘too far oﬀ’. Remarkably, the

inequality does not depend on the complexity of the curves ‘visited’. The missing proofs in

this section are deferred to Appendix [B.](#br27)[ ](#br27)Lemma [22](#br15)[ ](#br15)is illustrated in Figure [4.](#br16)

▶ Lemma 22 (Iterated triangle inequality). Let s ∈ X<sup>d</sup>, t ∈ X<sup>d</sup> and X = (x , . . . , x ) be an

1

r

ℓ

ℓ′

arbitrary ordered set of curves in X . Then

d

m



!

X

dtw (s, t) ≤ (ℓ + ℓ )<sup>1</sup> dtw (s, x ) + dtw (x , x ) + dtw (x , t)

.

′

/p

p

p

1

p

i

i+1

p

r

i<r

Proof. To ease exposition assume, that r = 2, that is X = (x, y). Let W be an optimal

sx

traversal of s and x realising dtw (s, x). Similarly deﬁne W and W . From this we now

p

xy

yt

construct a traversal of s and t endowed with additional information on which vertices of x

and y where used to match the vertices of s and t. More precisely we will construct an ordered

set W of indices ((α , β , γ , δ ), (α , β , γ , δ ), . . .), such that for any i ≥ 2 it holds that

1

1

1

1

2

2

2

2

(α , δ ) − (α , δ ) ∈ {(0, 1), (1, 0), (1, 1)}, and for any i ≥ 1 it holds that (α , β ) ∈ W

,

sx

i

i

i−1 i−1

i

i

(β , γ ) ∈ W and (γ , δ ) ∈ W . Refer to Figure [4](#br16)[ ](#br16)for a schematic view of the constructed

i

i

xy

i

i

yt

set W.

We begin with W = ((1, 1, 1, 1)), which clearly has the stated properties as (1, 1) is in any

traversal. Now recursively deﬁne the next element in W based on the last element (α, β, γ, δ)

of W.

If α = ℓ and δ = ℓ we stop adding elements to W. Otherwise, if δ < ℓ let δ = δ + 1.

′

′

′

From this let γ = min{j ≥ γ | (j, δ ) ∈ W }, which exists, because W itself is a traversal.

′

′

yt

yt

Similarly from γ deﬁne β and from β deﬁne α . If α ≤ α + 1, then (α , β , γ , δ ) is added

′

′

′

′

′

′

′

′

′

to W, and the steps are recursively repeated. Observe that α ≥ α, which implies that the

′

properties of W are preserved. If instead α > α + 1, let α = α + 1. From this deﬁne

′

′′

= min{b ≥ β | (α , b) ∈ W }. Clearly β < β , as α < α . Similarly from this deﬁne γ

<sub>β</sub>′′

′′

′′

′

′′

′

′′

sx

for which it holds that γ < γ and from this deﬁne δ for which it holds that δ < δ = δ +1.

′′

′

′′

′′

′

But by deﬁnition δ ≤ δ , thus δ = δ . Thus we add (α , β , γ , δ ) preserving the properties

′′

′′

′′ ′′ ′′ ′′

of W. From here recursively repeat the steps above.

In the case where δ = ℓ , we set α = α + 1. From this we similarly deﬁne β = min{b ≥

′

′

′

β | (α , b

) ∈ W }, from which we similarly deﬁne γ , from which we deﬁne δ . But as

′

′

′

sx

ℓ

= δ ≤ δ ≤ ℓ it follows that δ = ℓ and thus adding (α , β , γ , δ ) to W also preservers its

′

′

′

′

′

′

′

′

′

properties. From here recursively repeat the steps above.

![ref5]![ref2]

<a name="br16"></a> 

16

Fast Approximations and Coresets for (k, ℓ)-Median under Dynamic Time Warping

t

y

x

s

<sup>t</sup>1 <sup>t</sup>2 <sup>t</sup>3 <sup>t</sup>4

t

W<sub>yt</sub>

y

x

s

W<sub>xy</sub>

W<sub>sx</sub>

s<sub>1</sub> s<sub>2</sub> s<sub>3</sub> s<sub>4</sub>

ꢀ

ꢁ

W = (α , β , γ , δ ) (α , β , γ , δ ) (α , β , γ , δ ) (α , β , γ , δ ) (α , β , γ , δ ) (α , β , γ , δ )

1

1

1

1

2

2

2

2

3

3

3

3

4

4

4

4

5

5

5

5

6

6

6

6

Figure 4 Illustration of how the optimal traversals W<sub>sx</sub>, W<sub>xy</sub> and W<sub>yt</sub> of visited curves can be

‘composed’ to yield a set W that induces a traversal Wf (in red) of s and t. Any single matched pair

of vertices in W , W or W is at most |W| ≤ ℓ + ℓ times a part of W.

′

sx xy

yt

Observe now that











1/p

1/p

X

X

∥s − x ∥

≤

|W| · ∥s − x ∥

\=

|W|<sup>1/p dtw(</sup>s, x .

)

p



<sup>p</sup>

<sup>p</sup>

α

β

2

α

β

2

(α,β,γ,δ)∈W

(α,β)∈W<sub>sx</sub>

Similarly for x and y, and y and t. Further, we aquire a traversal Wf = ((α , δ ), . . .) of s

1

1

and t from W by ignoring the middle two indices of each element of W. Now overall





1/p





1/p

X

X

dtw (s, t) ≤ 

∥s − t ∥<sup>p</sup> = 

∥s − t ∥<sup>p</sup>

p



α

δ

<sub>2</sub>

α

δ

2

(α,β,γ,δ)∈W

<sup>(</sup><sub>α,δ</sub><sup>)</sup><sub>∈W</sub>e





1/p

X

≤ 

(∥s − x ∥<sub>2</sub> + ∥x − y ∥<sub>2</sub> + ∥y − t ∥<sub>2</sub>)

<sup>p</sup>

α

β

β

γ

γ

δ

(α,β,γ,δ)∈W

≤ |W|

<sup>1</sup><sub>/p</sub> dtw( ) +

s, x |W|<sup>1/p dtw(</sup>x, y |W|<sup>1/p dtw(</sup>y, t ,

) +

p

)

p

p

which concludes the proof, as |W| = |Wf| ≤ ℓ + ℓ .

′

Observe that this extends to r > 2 in a straightforward way. In this case W consists of

tuples of length 2 + r, while |W| = |Wf| ≤ ℓ + ℓ still holds.

◀

′

▶ Deﬁnition 23 (metric closure). Let (X, ϕ) be a ﬁnite set endowed with a distance function

ϕ : X × X → R. The metric closure ϕ of ϕ is the function

X

ϕ : X × X → R, (s, t) →

min

ϕ(τ , τ <sub>+1</sub>).

i

i

r≥2,{τ1,...,τ }⊂X

r

s=τ1,t=τ

i<r

r

![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCABvADsDASIAAhEBAxEB/8QAHAAAAwEAAwEBAAAAAAAAAAAACAkKBwQFBgEC/8QAOhAAAAYCAQIFAQUFBwUAAAAAAQIDBAUGBwgRCRIAExQVIUEWFyMxUSIkNna1CiU5YXF4szIzNIKR/8QAGgEBAQADAQEAAAAAAAAAAAAAAAIBBQYDBP/EACkRAAICAQMEAgEEAwAAAAAAAAECAxEABBIhBQYTMSJBcRRRYZEWMvD/2gAMAwEAAhEDEQA/AHF9XLq6560m2Dq+B8L12KsD62w9Rn3c/bayuet0iEWePG8s7WkE3piWtWdBZsK3JYUlP9rTFQ8t74UGNFdelXMpEwbh2cnqHsO1eOQRIKYC6UboqnUQ5MYEigZQ3Acn+BAO79Us9Q7A9L2T351BwHdUXo1DKmne/lYl0ot77O9OzWca4ikLeVIi5OxcNF/L9OsDdYSEVX7QKI/JzaGZvtmWMHNqrl40UOw+D7LM4Mz83rsaRhWlsl0EjFCdlaS2UWBzIY/fpvmoQU84bsPczIPRSZEFA/jW6HR62HXzySSrKjhdi0FEIWNdwFfKQsxLlnLFTagBRtHV9ydx9t9Z7a7Q6P0ftE9udc7daVu5O4Wcv/krM2nESoPBD40hSCXbGskoZ9Qxdi3Lhr1h+oLtnoHUsNZE121WldgMbJWyVs+zN4a+cnBYuxJUVoJOVi5ZZKNfFh39vJOqKR1qVOonEFrrvmJkfVfuzm4yXLIsGEmdIjYXsfHSZkF3HcRBR02IumRAQSLyUfNMUHHxz5YD5Qfl48llPGlOzJjm34ryBWYm20i+1+RgLRWrIyRk4uTjJBPsFq+iVv3eRRKcCq+nUVTKJky8mDnkAv6eOQrc3ql11MzJa5e0Z21KnEKNP2C1ySy92yzi5yLtHFmfJxkcrhOKjsqhE2dKLiAlpc7H7LuQO/WBQvb96eS5PKQT5GqufjxQuz6zlC0bs8kQISRt4BqwSqg+gByRfoe759n51M90750/dVrds1QNe5zYtzULTWm1tpkVaT1hjXqjJhJBO5FkZwlesHp42tnbsizImjRL2v2xhcJ9vBjhwxYb3bMWUWz5Np8Lj++z9djZa00mu2Y1zhKxKvm5F14iNthoiANYWzUDlAkoMJFi45EfRpcfIBdaFJMnSo3tApADjAFsKAcf9JQUZcFD9AD6B9PDJaYQidTrSaZSkISBhylKUOAAAjm3Hx4vGDrvRsa41D1Fz/s01pp8gr4TxzN31Klpzg1taxqQ6ZDJxaM4EbMe2qulFCpkc+2POwTc+Sb6c/WXMGRMza/4jyrk/Dr7BuQL9SIaz2zENhllH83j6Yk0PPcVqTeGi40zl0xASAosLBoJ+8B8gn5eBi60pzE6VG9hyGEpi4AtQlMHwICCjH5AfDFKm3A9WrRxVXATQEOIgVQQABGOb/kHHwHhjFa7CILK9XLp3MhMPKuse9/mHSVM2BQxHmtoAftIU/YJe4QKHI88/I+OfnOINqVu7i3aKOFNri3ahtTtV9hTlTCVl2t7hnEgrq7Yo1s99BE1DHtaLPZdDKFgNMgucZKq+XFyHkmFv+c/GKHWA6dQiYAAusW+YmERAAKHrtavkRH4AP8AXwdWzGDajsvg/KGBb+1eyFJytTpWn2FlGyrmGeLM5BMo9qUk0KZdn+MmiY7hIDLJkKYE0z9wh489h8hlEkm4iq3fH6qgRxVD+h+2ZZneJYZHeSNf9Vdrr8Hg+7P5JP3mt95FECqoqIimkUgpuQ/HIkU4cnVIQAAVSlACiXgQA4D8G4D5V/uPHl1ZzniTfmBcKRVKjF2OG9vo2PMEJWH+E5xdNSv57yNJNDuHsqGtZmc+WpV40NIeZ95056WSivxgebJ0/wDMVzyfhp1TMxPGK+wevlyl8E5vcQkc3gK9YLzQyNEHdxpMUAoSCuNbCm7IWryb+LijSZmEkJWKYoH8FBkTG9PyzSrZjy9VmBudUusDI1ewwdjhWMvCz0Q/S7BbyrR6RVJ20bHAFVWyyKqahxTHgezkMogQMNzNuYuSxs2avn9uOB9Zj8ACgBQ49cYBfWZetnvSg3ldNlPMRea92lw3NwICoiqLBRM//sQ5Tfn+Rg8MtppgNVa6IDzxCRACH1AwRrXko/5h9eOQ/wA/E2uyl9t0P0i+pVqDl60z9rzHqTiGx0w1htsk5e3jKuFJNYEcQZ8swmWfN2ZsnhCWtJGJGXknceWtGBycAWS5pIpYgNTrgl7+32OIABOIiIgEa1DkvI8gUfoA8fUePnxeMXN1qP8ACm3t/wBv9r/5GPhj1Q/hOsfy/Df05v4XD1qP8Kbe3/b/AGv/AJGPhj1Q/hOsfy/Df05v4YyQzrEdUPJWgnWN00uFl11WyBgOgYayjAyVnqM6wmblKQmaJTHqWU35KqMikaFPjgtLqp4l3aRg69M/aR4BZQwsFAJXRFyRZSLYyBhMiR2zav0hUIUigIvECLk81FETkB0QpgKYCiZMDcgkocBMPhPuy+G8R2XrBaPGsuMMe2MclakbvVvIoz9MrkyN+r8E/wBeTQkJdRko1yNqioc0pJGio+d9e0jzSD0zRFEXa/mOUcNGjVqLdEhG6DdmCYGSMmkoig3J+CQDmMQEyJJgcEzGMUqJCj2iUOeWMm76nG+NL6TfUB10zm+oeS7zBbw4yuGLsr1PGMTF2y72CzYLfVUMHOa3CzctFMIdpBBlXIA2MkE9M7soSEb65q79pbeVRnGu0lkGb8CukU3rQJIqaphH0aL1FFXyFkTm7kz8CBQSQKoRuchwL2gbkVUYQolS3e3QvW22Ra1WbZjbUuyz2CtNo6wQrB69rV5j3DZXOmb66/Zs3latFXyIDTGYYptpZaSma8as2j2hOACUdjIt79MiUo/sAYQKYAOf9s4FNwIlA5uTgUeA5KA8D+nhjJZv7S+xJr7qlbNraguLSRyXXktOMp1hn57CPyNVMvLouazfb08jCHfWWVwaFVnvu4jpVB60afb61AmrHmdHB1STgrIFUynh7G+Q6NJrTVRtlPgZeAll4qVg1X8cvHNwRcmiZ1jGTDETiU34EjHs3RQDk6BeQEVm5lgq7vRu1F622uvwFv1g1BTquYM01m2wTSbquUc8WlzIkwxDQciVnJM5ZjiFOrXsct48szhhDOButNM/hZnsSFi2pErhs2aItPIZopEIHpRTRICaaZQKk2SSRKKKaPYHaYjbjySlIVAvAjwxi1eto49J0nt8nXlqL+n17uC4tkewVnQIFarC3RA5iAKqoJiVMAMAiYQAB+efBR6nZlitjdasIZzhajkCmxWT8cVm1sKtdECQtrg275gmBI+diiPlyMJBLyh85uVZQCdxf2h58I0/tCl+tOb8BZp1exxYJuKquEcLO9oNmrXRpmRiJ+pe0u2iWD8VWVgs4imlkpefA+8cZVszGbJHjj1EJluw9Wx9VR/U0SKVatHMucpjQEOIlKIiUP7vb/kIAIf/AAeP08MYsnYA3Z1e+nQsYAAptZN8yAAnIAir63WoQSARMBROb57QAR54Hxq/UJzBeqZjKMxdg2bJF7F7H2RlhLDz4zCPmU6JMWQp0ZHK9prDlu+fSuM6Un6YlzmWsPLIQRpiI8wiQv0hPjGwKpEOr107VHCYIkT1m36FRZcRFsY6LvWkfUF5EQbpJlE5hWMCZCByYTAADx2OsUattPtnlncmWAH+OsLKWbWvU0jo51UWpVnDdXPeU61OwpjQd/x3msGGMgqc0q/sDNiNIlPZVYz1r/1jGHjr5hii654jxphTG0WMNQcY02MqVWjwfSMmRrGRyfHlDISSzmQeGBZRU4KSThRUCnKCZvkwB0O2mxsFq3gDI+Z5aHmbS7qkG7NWqPWAYrXW+2lZMxIeo0iIkF25bBaJNbuUjoZr5zl6VsuCKKvlmDxuhPOIRUomIJUE+CkBM4dxz894rFQKBR7O0OE24d5QEe4oAIeFZ2szzbfqAVfGzRRd9gPRo0Pla9rqESdV287ITjlwGImFRtcSVRZvaNf063dPvMp8jKNAMGRaqMvEvuxsLdjCd0Z13seAcDRcdkmZiLVnDJU7K5g2DukIm+QgrtmW9JsVrlZIeGkWrFOusJIGEcmELGxcTHNDtTA3j0O4wn2POmYKFr5irIWZ8mTyNcoWLKs9ttqmnTCSet4uNZJmEjoyMUzeSLkhlhTTWSj0HDkAOUQT4+Q2kFUwD5UIHaUDD+2X4KP5CPz8AP6j8eFSbPSR9qNxMSabRBEn+NcHOals1tOuUFkzIL+e8SwBimxwk35cBkDGmaQZZONeI6OY2H2UaTC++DDe5R3r2MBXZHFFzqPRY6gmZ841sIPYrZbGl9znl6uuX8fMvceubMaNWhcMQ9uYvJCUksb0BMHylMjHUu6aRYz0v6NJAztx30PU1sCVRq6aZzFISvQwFKPKggHtzf4E6gmOb/UwiI/UfCzusygkl0rt6CJlKUUde7gh+ECKaYETPGlSSFNPtExUih2pE7RFEvIGAgGLyzqofwnWP5fhv6c38MZI917OoLV9R95NTIJtaL1iy4PdedoMUNs2K4ytM9Tcay2wD7DDauX1tILV2UhMjJ1olJlS2umUkJ2zQhZGJPIxTYJJiK9RmAMX4/wRhbF2I8Vxrav47x7TIqvVaHScSblNGIbtU1EAZLzDl3JOExMsdYou3C6xfN7TmDgChJN1/wDpTSm8G+2Fpum32/yOU7vqPsJbKBix1Z67GY0YXXXB7jB5VzMzTTP08ejcgyDI/b6QSforvE4KAKxdMDJHFd/nT53wr+zOPseQ14aRtTylNY2ir1T1mUTP1yl5pxu4QKyDJuHmVwcO5ZRg1XZgW60R1JTVtxoD2uHt/tpbRCi7Yz1PUg6jmAemxg57lbNllKpNy51q3jShxDhua/X63SBSpINoWKERSLHtV1Gp5GfeswhYjzWxZF439YgVXRtIcBTOvuvUJE3eej7pmS+Sj7K+w2RYJo/ascmZguCLNxcrozhDpoowp5gjaNIeNho2LYJmbB6RikAn7un2+6eemm9H2DcbbYPrOWZLFPvK9Dfy1luEEalL2L200pJoFrlng0H3qVYWKVRSmk5AhDsOECFKo4Kp5PEGWMia73Ok6s7LTEhZVZ4EonXvZSVbMo+Ly8m2Mm3JRMoqwzOOq9dzmiktHizjWDWsx+TyOHY46rIjVbGJWMIzZfZjDuqODLjnzNd5h6DQKjFecpN2BdBivISqiLgzCvsYxwKTuUm5BREwM4OObrSzsibg7ZqcEVTFHPpuUaTa4RWz/eXVbmMv7W2WQz1kh9WLHC5BqkA9vBUFovHWP8iQb2bbT2KaURs7NUm6VilmcaaZlARc8uVAHe9sdScA7q4nUwrshjOJyjRQsURb21dmJKfioxpb4dGQJAz4LV+XhnzxWHVeLqAxM8UbuCrCRdssAlAB56b+Qp77m19acjvUHOcdTp17ha+uCRELSlrxCQIERrOaoLH0fHRCdUxvkkicknSBPEps5Ua5NDFvHZmbsUmM8z1oHCbDpRb4ySx0Ekmev11fSCqQGUSMCfo1XKjbzPMERAiQiVPkQMIAAFEePBZ6v5Up2xevGGc44Zuc1J4sydj2t2ujP3sM4iXTqAfsEgarLxtgjkZpmoYUz8oySRHJeP2y8CHgTOtU2bvulHvazcoEM2c66XRFduYx+06LkWKaqRFEDkTMUSCYpzAInADAJDFH58HPrxhbGODcGYnxFimqoU7HOP6JXq7T6uxkZh40g4VnHpC2YN3UrIv5FdNIVDiCjx45XHuHuVN8cMYvnP0ekbq79PBFQxnBFtXt8iKJueBTUMV5raHeYqQJfBymMCpCiCZg4DsAA+cUd4Jq1U2nyTppZJCTokVnKQtm4WjWWagsxPlbFmVIlZinsRFVlcrJ3UahjXHZpXEZsf0SxVB7B2L3uxBKoWU0W2GL2/PCgp9XXp1rqDyUNZN8S/P0E73Wzj6fkPaP+vH6h43XqCYYueUsPN7xhaPYOtitep1hm7BBZSRbQsBLX+nlXXj6rdpESIvnFBnE1VTWWFaSsSSTMxYAu7KKCXFGGSMnc28A0aoett88kEgix/I/Jj4M+q00eoDzwgs8m0AQC6oiqaqI9e7+85WvmxdhfXJbXLZGPhq1sfDQ5paOd15m9aY9z3S2JgQc5KxESVeSrpw0bHM3PeaV75YJzGp5OvDZHbYtmiiqkfljFmPs1USYxdk2ttZ+o2BIrVyxM5et3TRyUihUHkbNRjppKwM2xFQVI6XjHzGQbGFUEHBSmUKYXIllh7qJ6x4ozBV5ixV1K216EyphfKdaZFq2Rcd2RygoEXa62jKtpMWLtN0it50BNM5WDfpIpBJxb0yLYyPYa755thrqrrFs0nCwGyEFX1JGGk4Vg8hKLsVTYQxUJPJ+KkJN/JLj7WZywG80RvLz0zjc0zXxsUsqWyxhQn3/AN7/AJxHIJEU7CpA2sx48hWrkA9BWPoDPL4ny3edfL1WNX9k5+QtsZaHqMZrvs9IFZoJ5SXAwINMe5VdMGrauQWbCJnaDDnas65FZWKrJmpFTbGqU6Y2e7Ppp6z7RYn3bryysRjm/FjNfdwCRyZUmTqtPVyqYczhkC0SJHsRV8f689+QivXSCcK2lRyMAyUmp6Fr5bAco4hx7mSjWHGuU4Jra6baG5kn0W8UWbOmoplOCbiOlYxVk/iXceKvewkY5w2fsznP2ORKcwGXjLoOIKIktBd3khydiDO6TrHmDs8W3ubNMospQgN4jDWYl4IYZpH5na8NV4eUYLV9nlYp5MKtUYsarKmcsvO56yrhF30pN53LUSOGamvdsM3cpiIJHKc0ecFUDCIg4buCGKZJYgmIJSiIGEDBwzaofwnWP5fhv6c38Tm7L5OuT/oz9RXXPN9oXsuxurOH7fjLIk7MlbQ8/lGttlGqFI2CZ08iSLqv0fLBWM2jUTKC8ZvBqsqDOSeC2X8uiOqSkcjV62kq+apqEgIcp0zrEKcpgjm/IGKIgID/AJD4YxJu2mwuE8cdWTSdW7ZJq1ZLRNc9tEroSVsEeyUq6uS32EAo60ydcxQZoWk1VsBYFRUqZX/tEj6cFPTqdrqW7hnJM03BSpuUn6IrNAFUqqS7BUhTAoU5SlIVFQggJgMCgEAQEREB+UgdSLo6VPqC5br+YJHJCOO7zQm9VLRZqGq6q6oM2rl44sEfkAgTTY10ZvjJRIVdJqpXwq/kzXqPefeEwYOxhkSRlfYgdwZRGKas0TqglwKjdqgVJRMqPcIlBTgDf9w36fTkdV0/Uax0B1GxhvUxsrk+QFFDblPC7So9Xdjnih2PdHTe09F210XV9CfXJ1mdCO7W1ap+nkBMTxHQ7I42TyCWdKnkZqhV+Fk5WtgR8pq/vRkjV9VciOJ9rWFp2g13SRIEtJBfYNyyT2ignbpmZpHVGiQBZ/EX3W1wsYY/Mhax9wfeXwib2eNfKTsLUyV2zLTUFN1+SbzdAyDUXaMXc8Z3GLKoMbZqpMqtHoMH6CqwGctHbd7HuSppleM1zERMkE/UY6U+sPUlmMMSOfp7L9dlMNOJ8uPXmK7ohUP3u7rQC8kSeA0JKLP2ZjVOO4RRXZhwCgHE/cUSNHjGiTJlHNUFTqIpMkmbMTF7Fxbtkk0SCur3CB1QKmUTH7A7jGEeA8bdjZHscD2bPPPP95yLhQI/GQYTEhgq7EVfEOaFvXLEcGx92AFWv2wV1b3Jvq3s/wC1w+yMVASs1X7DEMF4SkbB1SvGaITGQcSxz91IuCrQgvY0b9RU5mxP8dGnK2aYm3n2lYlSJjJuJ6LmOgT2P8oQLW1U+zNTspKIeN1kkkSGKcqblBVNQr2LeoCYThINXSLlEwJikomHeBxw3u0Ywl1AsPkwTneQyJC1aKt9avqFhxbaG9OuzaVggkQjWqFgViJg5IxczpU8k0TbEK8Mi1FQweSUBKTEdBQxxjSk0T7V3LIAVOBjoYl1yLLJT15spY5MUkZa0zDdlGt5GYWJ/wCS6SYNSKCACCJfy8TkZGp1s3OdNEMX5DsFyQueTcfZ6xi+1ONlmBbmWk8qV60PWH3HG2cupGjyEbX3XJJtdvsm+bxFfbZnJkuxDDRFWGlPRkrM6mwcKVatKEclalPAQ5itzpAcyIDHtx8sxx4Ewl+oiAc/p4zra/XbGG1uu2W9eMwpS441yxT5Ko25WtuWsfZGsTJJ+Uu6r8k7YySMZLoAPLN+LFwZuImEqY9wh47PWnAdN1bwHifXfHDmwyNEw/S4ij1V9cJYs5Z3MNDpCk0VnJcrVmWQkDEMPnOQaoAoIBwmXjwxn//Z)![](data:image/jpeg;base64,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)![](data:image/jpeg;base64,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)![](data:image/jpeg;base64,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)![](data:image/jpeg;base64,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)![](data:image/jpeg;base64,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)![](data:image/jpeg;base64,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)![ref3]![ref10]![ref10]

<a name="br17"></a> 

Conradi, Kolbe, Psarros and Rohde

17

2

3

2

3

2

3

2

2

2

8

8

3

2

2

2

1

1

1

5

5

3

3

3

3

1

1

1

1

1

x

1

Figure 5 Illustration of the metric closure. On the left a distance function on ﬁve points

represented as a graph. In the middle the shortest path tree rooted at x inducing all values of the

metric closure of the distance function from some element to x. On the right the metric closure.

The metric closure of any distance function is a semimetric and can be extended to a metric

by removing duplicates or small (symbolic) perturbations. Note that the metric closure of

dtw can be strictly smaller than dtw because DTW may violate the triangle inequality

p

p

(see Figure [3).](#br9)

▶ Observation 24. Let X be a ﬁnite set with distance function ϕ. Let Y ⊂ X. Then for

any σ, τ ∈ Y it holds that ϕ(σ, τ) ≤ ϕ| (σ, τ) ≤ ϕ(σ, τ).

Y

By Lemma [22](#br15)[ ](#br15)and Observation [24,](#br17)[ ](#br17)dtw on any ﬁnite set of curves in X is approximated

d

m

p

by its metric closure, with approximation constant depending on m.

▶ Lemma 25. For any set of curves X and two curves σ, τ ∈ X of complexity at most m it

holds that dtw (σ, τ) ≤ (2m)<sup>1</sup> dtw | (σ, τ) ≤ (2m)<sup>1</sup> dtw (σ, τ).

/p

/p

p

p

X

p

▶ Lemma 26. Let X ⊂ X<sup>d</sup> be a set of n curves and k and ℓ be given. Let X∗ = {τ∗ |

m

τ ∈ X}, where τ

ꢀ

is a (1 + ε)-approximate ℓ-simpliﬁcation of τ. Let C ⊂ X be an (α, β)-

∗

∗

approximation of the k-median problem of X in the metric space (X , dtw | ). Then C is

∗

∗

p

X∗

ꢁ

a (4mℓ)<sup>1</sup><sub>/p</sub> ((4 + 2ε)α + 1 + ε) , β -approximation of the (k, ℓ)-median problem on X.

Proof. For any curve τ ∈ X let c(τ ) be the closest element among C under the metric

∗

∗

∗

dtw | and for any curve τ let c(τ) be the closest element among C under dtw , and let

p

X

p

∗

P

dtw (τ, c(τ)) be the cost of C. Let C<sup>opt</sup> = {c , . . . , c<sub>k</sub> } be an optimal solution

opt

opt

∆ =

p

1

τ∈X

to the (k, ℓ)-median problem on X with cost ∆ .

∗

Let V = {τ ∈ X | ∀j : dtw (τ, c<sup>opt</sup>) ≤ dtw (τ, c<sup>opt</sup>)} be the Voronoi cell of c<sup>opt</sup>, which

i

p

i

p

j

i

we assume partitions X by breaking ties arbitrarily. We can sassume that all sets V are

i

nonempty by removing the elements of C<sup>opt</sup> with empty V . For any i, ﬁx the closest

i

P

P

π ∈ V to c

<sup>opt</sup> under dtw . Letting ∆ =

dtw( <sup>opt</sup> ), we have ∆ =

c

, σ

∆ . Let

∗

i

∗

i

∗

i

i

i

∗

p

σ∈V<sub>i</sub>

i

i≤k

X

= X ∪ X ∪ C<sup>opt</sup>. By Observation [24,](#br17)[ ](#br17)for any i ≤ k and τ ∈ V , it holds that

i

dtw | (c<sup>opt</sup>, π ) ≤ dtw | (c<sup>opt</sup>, π ) + dtw | (π , π )

∗

i

∗

i

p

X

p

X

i

p

X

i

i

i

≤ dtw<sub>p</sub>(c , π<sub>i</sub>

opt ) + dtw (

)

∗

p π , π

i

i

i

≤ dtw<sub>p</sub>(c , π<sub>i</sub>

opt ) + (1 + ) dtw ( <sup>opt</sup>)

p π , c

i

ε

i

i

≤ (2 + ε) dtwp(c<sub>i</sub> , τ ,

opt

)

![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCABvAhIDASIAAhEBAxEB/8QAHwAAAgEFAQEBAQAAAAAAAAAAAAkHBAUGCAoBAwIL/8QASRAAAAYCAQMCBAMEBwYEAwkAAQIDBAUGBwgRAAkSEyEUFSIxFkFRGSMy1gpCUlhhcZUXJFWRlrEYJjOBJSfBNDU2Q1NiZnKh/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AO/jo6OjoDo6OvDDwAj+gCP/APnQe9HVs+O4Ic5jgT0yeqqmqT01EEz8+n6oAJwDjxNyPvyH5AIdU7qYaM2C8m5eIosEUvXO8WFJBsmh7iKp1l1UkiplD3MdQ6YAHA8h0Fe9EATTATkKUyyZTFOTz9Yo+XKJfqL4mPwHBvfjgfYerM4kUkUiKrizQSMCxnvxLsqRGgJeHrnEVClIqKfmXy8zogXkPqHnkF/5I7hlPlrXcsLaqUiwbY5zqT9as2Gv0cRjsYYwuaonThI3O+UDtpBPHtemlEJECWGLr1uKiSLen+CU9MoHj6K1Iz/tUUs7vflc6NBmQRkmmnmE3TiGxcxhZzz/ABNjjO94VVXHZmDSKzjEI2VLU8Xi3IaUMMcf5oBWoZJlnuGUhtL5CxDq5VLLtLnSlqTtftFZx6QXePsMW4qImgmOfMiopuy4ogpZVtJCykEa/Z/WJESYAgHw31po7Sveg2/ydWMiX3uT1qBgKC8aUu3Y/suF8XOW0JhLFUqtaChmDZCRXsSpqti28FYpBjq1Jt3hZ/8AC1xODNr8vAqrbd8I1nrRo+7wBquyZ4QuuYH9Q1S1XNS0i1uKp2TcgFlEaIV9LIlcOK5XmCEVLJup70pVVuVcoFYn8hHqOWuJafo3u9rUXHtQiMea87NYRR1syios0RVocLkrDho0dVMYUEgFTTr8neFMg5mUXZKoyBLMEK3MB435WIOgajjbLWOMzUaByJjC4Ql5o1oYITcVZK6uk+ipyJel80ZBIpVCqJoPQLy2UVKUTgRQAKPiPUjF8hMcHIKGKCRAcicnDZYg88AgXy9xJ4/WHIgXyL7j79LuyX29KQFqtWZ9XLrMap7ATsk/sU1bselE2Osl29AxFIF7sDjJFxGFytB11ReRMzgErLVgIWXkyA/D4vlPC2G22w+s4DG734idK0Svh8rbbfYTaKzOOJCFZcfiXKWbsfmIi61zgAKeNXi4RpZspqLENKlGRL8uAXQNfDjgOPYOA4/y696iTG2bcZ5mpVfyHiW61+91C0x7OSrkxCvSLNZVjJFMaPelAwA7Zt3QJKiiLxoiocCH5SKJB6kVN8qfgxQ80xAAIcpA9Jf1P/TUKp5fQQniYFOSm48ij+odBdujqibqrKiXzEhPApgUIUAOBzDx4GIpyHBQ4NyHj78h7+3vW9AdHR1QunJkFCiI+KQJKnOYxA9MDgJAIB1eRFMR5HxL6ZvP3HkPH3Cu6OrItJLt2vruEwREhTHcm5BQjcpRKAiTgCiuAibgg/u/MOR4Lxx19zunRDCQSF8RMb96X6xTIbgUBMkPjyBigcVDCcvgIFAAMBuQC6dHWAOcjVJi7XYPbVWWb9s5VbOGcnORUa4TBsPi4cAkq8Op4AYxPEihUzDyPPHAh1SKZVoJEklS3ylmTUSMqVQ1lh0yqIl8QOuX/ezCBCCcgAUfc/nyHHHQSV0dYrG2VpMNVH8W+YyEcRVZueQaOWzpqQ7cSgdQirRVwmqCgnDwJ6hDFAo+Qe4dU7y5wkWdFKZmYqIWcHXSRTkHzNoY525kin4B24biZQoqlA6SQK+BjAAn5EvkGZdHVp+PWP6gokIqmVNJZMxT/vFQVA4gn6fjwXkC8lU8jAYAH6Q9h6r26oLFOIGKfxVOnyX8vHj6R/8A3BzwPQffo6OvDDwUw+/sAj7e4+wfkH69B70dWv44AIZU6gEBMhVFk1U/TUblW+pIFg5MBRKUDAbjnkf06+DuZZsWC0m6eoN2CJCrHeLikg2TRN9lTrLqpJET9wATqHIACIc/foK96YCppgJiABlilEhyefrh4nH0S+4eJzCAGA/1cAUfpHnkLK4k0kU01VzM26JyKqPTOHpECMxTEgLn8lSFIoKahygYTnRAoiHsIj0vzIvcLrU9Z7rh/UjHdj20zLTpd1UrGzp5hh8O4wviCyqbKu50yydrJlx+wlkmUyoymIyr25JRSHcJ+hwYhxwSK08zbtSQLFvpld29ps0VGSZ6e4berwWGI2u2Ehns/jbOlhVUdjs2hDOm8Q1hbcWv4u8StpByMCX5uCTIKrOHcXrjRzlHGWpVNsG0edcdtrHDz5qOiMlifEF2hGr1UsLsFkVEi4YxCU+AknMY7RgrIWWGDkUASbAT1ekVdjn+kJZ42tqOwtn35qkgxpONZ3GwN8zYcw46hcJ4dqc7F3FaVmM321xaXpqtHyjmNiQgZBNg+I8TQkBKmUCBy8TfKOZa9aaQ+umsbJlhm3ZosNE1J1na1BEldgK9Z7O3kZiHjZqTQBdeCrTWkUa1RZpcUpNX/e0Wvwvi7OomuCxVml9uLvi6u0yqzc9UtUd9dYlNdahrXQaURviyrZT17NSYnELycAsmLF4yJVLDfUYx+sybrwCZlWibeULJKrsw6KMZ5UoOZKJXMk42t9fuVFtsPHWSCmq8qnJR0vBTjcHkW+aOCKkOZs8REiyKyyCJ1ScG9IBAQCSCgcTqA59QxQImDr1C8N1hHkQFsHkPIEEo+QB/ByUOR56XLkDt8VVpZLVmTU6+2fUrOdqlH11skvS0hXxNlPILtwZZtZdh8Rou4kMonhiPJlGNjkrXWStzTcicHIgoBesVY7gZ21oA0ZvriKSYVCAEIhtt3hpkrP4ik67G/uJ7KmbaoJGrvXFpKugiV4OnNJnKhikkJBsM8b5UCz0GvdHUVY8zNjjMFMr2Q8VXKAvFPtkbHS1bmYZ4VZpKxsw3F5FSAewOmjd60IdduV42QXUKA+SRBIYAz9N6qp4nL9aQgBSGKQATceqAmSXKp5ckTKBDeYeJuBUIHP36C7dHVG3VWVEonEhPApiqJlADgc4iHicint9IABuQAvA8gPPt71nQHR0dWl5IGQUXRKZEihE01E/IwiIkEwAodUoF+ghBMUpRATeQnD2Dj3C7dHVgWmPSM6KYSlK19MypwAVPTSWMAJKHLwX6lAEBIQB+oPIwmDw8TVajl2Q4pmIUC+ZwFUnCgpgJgFDyS5IAlMmBhOfzDwNwHiby5ALp0dR84yRUGblZk7tlYZvWzhdu6ZyU7FR7lArQxk3DgUVHh1PH1PDwTUBM4lH34EBDqmUyrQSkSULfKWZNRIFiqHskOQqiPJCi4Jw7NwXyOQCkH3P588gAcdBJXR1isbZGsu0NIRj1k/jgVXQGQauG7lsQ7Y4JHMRRqqukt6xxASEA5TFDny4HqndXODjVUW8xNRUU4cmXTQSfPmTQ6h0FCJmACu3CBzqeRylMkiVUCmMBfUEePIMy6OrSL1cwqCiQihCFROQSm5UVBYBECCn4h6fAfUB/I3kAD9Ic9V7dUFiGMByn8VFCck+wCQ3AgPsHBg+wh7+/5j0H36Ojrww+JTDwI8FEeADkR4DngA/Mf0D9eg96OrZ8cBUzKqKAQqZAUWIoT01ECqgBkgVABP4iBQMBgDnkf04HimeTLRgwWk3bxBswRIRU7xwZFu2TROIAVU666qSRUx8ih5qHIHJg9uR46CvemApUuTEADK+IpnJ5euHgYQRKPkXwOYwAYD/VwBRDxHn2sriTRQSIq5Fkg3Mmoo9Fy8IkRr4nAi3PqEIRXhUxUzCZREAE32ER4Fe+QO4dW7NY7libUHHdj22zBT5d5U7AlUlDQWFcaXxm4VKnW845cUZyv4BbSDZjMKsJSMqlsRXcxZkASKCpVQwuK00zHtIQLFv1ld5YalMglIN9QsSuF67guOrk6QZORx5m1+qu/HZFWtyhI1KFvBYnGfqmjjOxriXxnotwo859yaqtoHMVM1AqFm2ly5jeJvcZZHlLbmd4Zwve6ug/ffh7OmT0k3IY/UBpDy7mLjUq7OFl1YwY4HLQFvi01b9pbvK7V5CxNbch90CssKfFSClCuNWv+KMUOYjHWCcU26uvZat27ZaUXsrw1LisqNHUJNYqfEQkC2KEc/GARsBvDponcBCPwtqpS9VtfWrbFVz2KvGPNXtdWdYQJX67CuWzZxb5SCm5VAFlYCtq4dx9eq4LwW8oo9M9axxkUivzrJY9E0inaZ9wfH9GplXh8d6w7b4EiMbLIyDFI1VDPWvzGsVHXrC+N+OEYUrjBzXJ02tXFUHpZVCqGkiPWfy4Gy4NBomR6PlKnQF5x9ZIq102zRjGeiJ2DVSfRczFTCCT2Nk0vTVAySUiksg8amVADHQER8QEREM5L5Cc4OAVMUATK49UvCC5vYQM3LyIiBTBwIc/SIgHI/cVv3zt31WEn7NlnUDIczqVmmclpSzzL2lInksK5IvDyQNzZM7YdReQoZPXj2L2WZxLdK3VkrBxI+uCqpUvSNjbHc7Net4GYb9YYk6lVYUflKO1+IwVt2FnlcjB+EkMm5mi/hIx7r8azSicapB0Vq8yYLX5yZiaxrgwF05Br3R1GVGy3QcpVOAvONrZB3GrWqNjpavS0S9TXayUbMMyyUZIcl8l2qLyPEXTdN2ggudMQA6aZgMUM3I+VMJTAHkQwAVMwEAE1vUL5prlU8gEqXgAgcPAeDmKHPQXbo6o26qyogJxKUClMVQgB5eR+fpORTkOS8APt4/mA8+3vWdAdHR0dAdHR0dAdfI6wEMJfE5hAomESgAgABxwA+4cCPPsH58D+nX05D9Q/wCYdYzMyLSKB5IvHDKOatEEnTuRdu0WTUqKHn9Mi8XEiLVuXyESnOYxB+sTeHAeQZD65REOAOYB59ygAgAh+Q+/36pH8owjWbh/IOkGLFsQ53Lx4sk2at0yh9Sy7hVQiaSRf6yhjABfz46W/kDuB1+RtdxxJqPj+x7YZ0r0o6rEmxrDo0ThLFVuMIEgiZ0yqs0evaHVZhROQ/8AMFIpeTVE04t2b5YoKaILYWhpXnbZxuzmd8Mw/H1aSMEo906wyDqBwhFtJbgbFjjL9oXcJrbXVISt2KEZI2Oi4vVTIWTEYgoyhyoBmWQ+4HWH1pt2GdTKZYtrtgK2/e1GehaOudhivE9zVMVOAb7CZTBB+vj2uS6iMkDOzwVVvB/CLlDfLP3BAVxlhphnjaZqlN785bCRpMsQrw2n+HHDmIwo0jpYObLjXNVkOqQdoqqqDeOQjZCcp2PhapklDBFm+anK338xfifH+H6TW8c4wpkDQaLU45hBVGnVtkjCw9cr0QU5WEdHMGyZk2xGhVjgVEpzBwYAA4CHvKrcxBTDxDwLyPAcAUBDgOBD39wH8jCACP5gHQRfjjDGOcQ0at43xdTK7j6iVGJb1+tVCsR6DCEgYFsHijEw7VAiZGTRMvAEIUDgUCgAB1JQoKcGEoiAj4pgUynmmCZAH954cF+s3P1hz78B7+3Vb5F/UP8AmHVHIPm8cweyDgRFuyaOXiwJgUyhkmyB11ATKJigc4kTN4FEwAYeAEQD36BUOQxX2K7nOJsSMDnksa6U44d5uzXTrMYSwSmVcvOmRNVrzTGIC5RmZ6jEx3l8qkm5LGL1Ys6gDD48JhwLWR+5pQLXdNVrHb8PVuXtec9eJ+A2IwJWohEHzaUy3jkHylRTma+c6CdmiDJScmdeAUeMSPwIImdoelwfDe2Imhkyj5b3SkESzLjcjKs3kXFNwlAAb+31cbqCOEMa3IB9VOFb47PJXIzCsMpKYjoMbC7Fk/cfGLAVlThsjMRTlg5I4cN3CLlisQClMoZB0kZAFyKiconAhDKAmudMpxEwmFMPH3DCsPZPp2bsWUXJmP7RC22qZCq7GwQFjgh9eIlWjlumJ3TI4+HrNwXFQhRHw4EOA+3vJIRwmRFJcpHBFQSTXI4EFzKpE8/IFTnLwsUfIOCGIAF4+489LB7Xq6OKq9mvSOSMEC+0+ynNU/FNIkVBcWptqfNnOXAN0sMkBSEnht6kJeis50ybRw/GHceswa+iT1Gq+RR9wEBD788h9v16Bc+TO3TQnFwtuadZLdP6j7D26Td2y033FBfhaZl28ibzh5jYvGzR1Dts0xMCopIGYQMpPwyZSzEqQr1P4w5gjNptfsXqm+Uhd8cVoS2L4ZM0cw3Cwg1VlKApXofgs/k3Ymgqt47/AMPEU/F3Fmh4KAsGWgMcJVMX5QZEO5bJyH6h/wAw6ssqkRwZJJVsk7bqeRFkXAAZqYo+P1KAJVAOoX39JMSeJ+T+RycfUGAYjyvjjMdSgb3iy4V27Uifim83XZytviLsXEBKJlWh3YoCUjpuR6gRQ6JXaLdYQTMApgIGAJQO9bpnBM5hKcwHMUoh7immJQOqHv8A+mAnKHP3+oPb79LnyB278fq3a05j1oull1HzvcJRW13C2YxYkPQMs3lZRRWMsGxGJm8nAMszsoQV5UGULKWeHImMxIiV6n65gGOW+12ymrMo8h96cSLz+LGCibdluRgMqthoa9WrnqN5vJexeL3rSCca6tJgz+LcxtTo9gzsUeZRE8uARiKr4GvC8QAAEwiQBP4F8g/i/Q4cCPJDf1R9hH9OsesTpyiwk3LGPO8fs497IRqAkKPxj6ObqKNkCGEeU1HKpgQbnADiPqmN4gJPE2AYry9jDMtIisg4qusBcqlMxzCXYSEAcq5nMI9TMpHfEMFTN38Uuul5nKhINkHBQIcPS5KYCyY4Aqjd2U6SbhEhfEU11BBuJSiIFAg+BxA4ciKnBBATeHAj9wBAbe+/0gHZecwHZqBiPVTQrF1oMrJZpa5FtstsLmOvQkuDRWKMagPaTjOPjrDAERfIvqsW2LN3Kr0ofOUvgwFeZsbdq3Zl1csq3LaDuqbl5XVv84hO1qtYOt8trJR8eEUO/UloGLq0ZZchN3sC5O4YlhmIOWBIBFiqil8YD0TIOcjy8HACHX9IhBKALKCdc6o8ev8AFAIceoQQIBDAY/PJvt1eOg5tsjf0Wftp5bu1syRkuf23vF5u8y4sdts1j2NlJGXss67UOq4lJZ4tXBUXdqKqqKHVP5GOY5hEQ/PDA/oknaeATpgjs+DcSpgBT5/kVOTpgYEziUa4AFFADnKQAEQEqhv4eOB6fOjoENVXsSY615wLa8OaW7o706snffNpeovK9n2bnaNWLlMnQM7sj/HaaEA1niq+gUXTA01Hg6MVIxnRPTADcjneQ7X3e+tGY8U4TjH2yvcix9gKHlrlj3ZlxSlGtmWnsquYV9ZqdLPVbnOrPWdScUyGGHAXgAzF458USer1/S9eeYJFMmUpzlUKYpTHEgCIAb8wKfnjnnx4AB/UPv1aipCBk1B9MwFDlJyCfLg66nuoUERHjxWEoCAmV9hKA8foHLh2bN3u8G1WxRq93Ge33sWjHtY+SqwblSifg4WcNvlaNQQytALgoJxSaJS6lgyCFidvDufhylr5wcGMl1MxrsijQihRUVKcDKeQF44D6eCiQR/dGHn2S5MAcD9Y9U4D6hRRUN5CdZMAarB6R/TMCggKwlFQFBUEAExfYoCX3Hq2PpBpEIvpF87asGzUCu3si9epM2REkvIBTePXBk0WzInn+5McxkxDy8xTHxAwZX65eQ4A5gHn3AAEOQ/Iff79UchKx8YycP5F0gwYtymFw9eLotWrcgD4+qs4VUKmknyIB5mMHHIe3S3r/wBwOFlrTcsTag44sm2Obq9LvarJp150aFwTiq5pqmTiWmc8sLNJCRotemARkxRm6RRsm/8A3U5/3ISgkZTDUNJc1bOIM5vfXLYz1YkR+aPNOsRFdQOBYdpLgCs7jnLkyu4Kfa6vsl0I5OBsNjpeL10RavXQwKZpM6LYMsv/AHAYSXs1uw/qBQbNthnSsykjSrClT1zROFcTX1u5FNhD7CZcK3k3eP2UmmymTxU5BU27g8+VPii0SDwOawMdJczbRtkZ/f7LS9kqssQj0dQsTOnMHr+wiZkAdz+Ms0OjqnJtDEMHiEWjAWubrGPlUE2D1wECQ0udJqwLGmLqJiWk1rHeNahB0Sj1GKiq/UKdXWiMPD1qsQDUWMPFRzFsmKTZGOaCm2TQKYSATgoGAC8DJ6BiCmHiAELyPiXgCh4h9hKHI8lH8h9hH8wDoI0x7h/H+KKVWcd4zqEBQqNT4RhWqzUa2wQYwtfrsYim3YwsM3RKmRkxaoIopIplKcCJpEIAe3PUimQUEDiAiAmMBPE5xOmCRQHgwE4AAOYOPIOeOffn26rPIv6h/wAw6tFgmmVdgZqfkPUMxg4mSmHhUClUXM0jGa71yVBM50yqLCigcE0zKEKc/BROUBEwAqa2x8xsZ3SKLUGbj5phfRfFD695epFnWMWKR2JzY4r8nrDkWjRQFct5eVptDqWcYhewLni3lZTtIMGKMglNulWsVd/Gj2WW7f05m3HmQrZjC4ab5Koe1tam6HGKO7TKK0AJquuqzCvUHzB1COJeKuz5dOTbg+OCTFRsZicjoyyE69rc5cj4iyHuM/AZdXdzJ85nPGtskv8A8dG1hmnL+V1wod6KPqkh1sYU6wv4ZpWGcjMxtaGScs46SeonMqLDLRBN7rT5+rPFzN05yIlob4t42I9TIlJRrqPB83ROsmDw7b4kFSnVFsZUR4EE/PkoR3rZnKk7UYGxPsJjUZr8CZpx3Wci0tSxxSMPNpV21xqEtHJzcI3fyKEZJfDrJeu1RknpCKlOALn9MBNOQR4mQFFcpHBFRRIuR0IOTKpJlOAguocvC/JjFECnKHHA+489ID/o/dkXxri/Z3t3WPI2RcqXTtv7CWjBprveoEK7GS2MH0nOMsR/g5oadnFUohGDpMn5w5VfgK8VRq0j5CYScKOEehIRAPuIB/mPQLiyP26aaFutuY9WL3ZdQc83GZkLnbrPi0BJjPL+Q3rg66Fr2KxM0eQjPM4xJXcySKYydhh/hfnkkYjkPXOUY4Z7c591WerQG+uJFHGPYZM8fHbg4SZKTeJ16pBCVnNZP2IqiyMSprqFgcLxDqEp9flMtpNwcSjT54csaRd22TkP1D/mHVjlkSODpIqtUnbdUDFWRclKZqIcBwc4GBQDrk+ySQkApyioJlCCUAMGC4lynjrL9Sgb5i23V66UefiGk3XZ2uPiuGLmvS6RHUE79ASkcNyv2YCsgV4i3XEhDcpAIGAJMO9bpnKmcwlOYDnKUQ9xSTMUp1fvx6YGMUOfuPkHt0ue+9u+ifjS1Zg1hvdo1EzhcZde3XGw40YkXxplm/O1lV2dp2Iw+3lK9H5mLDlczBI+Mk7PDlbnm3yhXgCqYpo7b7a7GatyT2D3txE7l8ZsFiIsdycCkVsuN3VUrYKMJbJexWOnrSAe66Gsbh9EO4qkUeYzqgmZeSansHjFIuH4NdF4gHHkIk8lBTL5AH1cfY4cCP7s39Uw8CIfl1jFslZKKhp2TiIpzMycXBzMpFw6ZvQGYkYxkqqyiSLgCogaVdAk1RH0jcGUBTxExPEcJxblrGWYqRE5BxZdIG4VGbjY2XYSEAcq4uYWQQFaKBwxVMg+inC6AicreQbN3KYEUKKPJDgWSnIFUbuiqJprpE4L4LKCCQlIbgopiBDiByc8mECiHqeHvx79BzqVyzf0j7Y/LGM7zF1nTvR3W+8QlYkbXR7eeSz3nTFaT6EUdSgy0U+ruNmUxYFJIzdF3VDyzBvDKA4STmXwIFMrsri7tU7IGseUrRs/3VN1cuvr5Z/xFWojDF5mNa6Tj9q6VfLydajKjG2HILd7CGXctAgm4PGJIBkwBikR2VwKqTnY8oAYPBRYUSkEhQWUE6x1fb1xdBxx6yZwAoGAx+eT+/V26Dmvvv8ARXe2XlG4WrIGQ5rbS5XO6zr2z2qw2DYuUfyljsMk4WdP5iVdq1wVHD1wu4XWWVPyZRRU5ueR6xMP6JJ2neTpgjs8CAgkJSnz/IKAJ0yiUigk/DoAUUSiZMpQEfIihvcv2Hp86OgQzXuxTRNf9f7HhbSvdrerVYyppGTpbuCz1N2Oh1OzTD9F3KzjvHBEq61nCvABcFmRpxh5rrJuRceSAEPyHd4Ltf8AfCsufsX4riVdme4nTNcIJ3N4x2lPTFGdnkZTJj2GtNgr0w8UuU2s+LS5avsGcIoZ9/uZEzFIkT1OA/plvfMCEFMoHOU/JQE4kDnxMHPIFNyJQER8R4AeP4g46tgJcGTUHwHx49FyUnLhRVQOVi+jyAACwgJ+RW9vEOffoOYvs5byd3d2vjLWbuQdv7YeJKzgZOtl3IkUAInLSccpHlq5MrV9fyFu4NDNpQZ+/BYJJ5ITJUCjBECTUVa9Psa7Io0TUL6ipTlFXzAvsAGEviUSCPKQ8D7JB5AUAEAObjnqnDyVL6ag+Xkun/uyv7sxSCU5iesIAoCgqCAHEo8B5F5E3PHVsfSLOGRfSEi9axrVqYjp9ISD1FkyKQoin6Tx45Mmi2bEFXhuCgiQQDgRIbxAQysFyiPAFMYB5DyAA8eQ+4c+X/04/P7dUUjKx8YycPpJ2hHsUCm+IevHCLRq3KI+AKLOVlCJIkExilA5zB7mD8+lsXfuCR1jsVwxbppjGy7ZZjgJuRqUk9iHh67rzi26snSiTeIzdl1VlKTFIjZNq0lzx8nS8e5JQdOGAICRNJUjgMUR0ay7soiynN9cxObnXn//AMTeah4rSc1bXqIYTRfjJnG2WQUeuC7TxcK/KyQr1zsNZxk6OWPO/NWWikidq1DKLz3BIawWG2Yk04x/ZNss3VuUkqVPHqi6kJg3FN7YvBIjB5/y6VrKvaCi/bMJk8TLwdLuxX6kcukKCBFQVCyMtH8wbPtkp3f/AC64uVWliEem1Hxco4ruvTKJmAB9M42zIUzlwls1HxEinGkr9xm69QliFjVnf4fRNJmSasGx1jWlYtptboGPKrDUql1OKioCp1OvtUomJrdYgGfy6GiY5i3TFJujGsQTaJoEMJSkAAKcALwMlNzFFMOA8S8j4gIAUPH8uAAR5Dj7D7c/oHQRzQsS0TF9NrdAxzVIKjUqpQkfWq5Va8xRZREFXYpBJtHQkQgiVMjJgzQQbpIJFKYE00SEKHAc9SAZBUQMP5nN4mKc/kmCZQEAEC8AACcOPMv25ER5Hj3rfIv6h/zDrHLhZoymVK0W+Y9c0VVa7N2STI0TIu7PHwcY5lHpWqB1USLOTNmqoIpHVSKooJSGUIAiYAV3Jxq2xfdEjogRLNYm0QxOEvdKpajctIXZrNJYGw4VyNQo7hyg9fQeHS5TrTuyKHjXsKlZF4ds2eN5RddHL+5xRJycwDDZlxjBy0/mbU/I9Sz1iOLYl+MiCzzEzqlWWQsdeMKKdiiY3Ft1vcmMed2zBuZqlL+sY8cVqvbO1vCSs3g627PW9VKwWbcnJNpz3V7s5XO5uUhrvbJaUntZ6jenCyYDHvcYYpn4qpNq83dyzCsgRSJjZF81SK4Owyy16PvVPnqnNJPHMTZYOarEmmkIJORj5yPdRjldB15HMRQWThYiDwUvMplCrCiAh4dB8KDfKflmjVW/48ssTaqZkCtRluqFlhjC4ipyuTTVvIRUswWEqfxDN6ydN12y4lIJ01CKeHP09Zf8uE6IpLkTcEVFMFyOfFwKhEiGAAWOcvDj6/AQA5A8eP8AAB6Wd2q5uZrGIb7qpcU04G2ab5NtGE65R3RACxQWuETMysdqtO2R0URJMvrph+Ci5xOf4QWmimVfuGLFRYWxGkAICHICAh+oCHHv9vf/AB6BbF87cFNi7Lacrai3yxac5rtc1JXKzSuMkfVxDlbIUs9UXdXLYHDTV9AsMwv27d9MoRhJKxxBmK0oo6IuYU/SUwNluPnPWF24g9/sNHh6JFEGPaba4YKtbMLnqcMYsW7yNn+McMYFzgJ9bpE8W5hqPAqZSbsFZVSM/EK6bIHi7Y+Q/UP+YdWKVRI4UTSVaou2xwEFUnJSnbDwACU4kMU4KrlMAAmkYpSGKJzCoUSgBgw/F+SKLlOrwl0x3aIC3VCciGU1BTVffJuWTiAmW6T+Ac+iHiu3B9GnScIkdJIL+nyBkim5AJDM9bkOCZzCU4lMfxEOBBMhgIKo+/AE8jFAB55+ovIBz0t65dueiRlts2VtUL5ZdPMx2uZdW60yWOI4j/EeUL/JrquFbhn/AAojL1uNzA8Zt3kujHN5G1Q4sHckL0jk5kPRVwJruBsJrJIPYffbDDprjpgt6bXb3BLlxeMSLVWBKeHXyFsFUn8ZVZTAEhcJFeKdRdApambGjJzJqMDWhVGPK9cg1/4xEPHzEyfkcSFA4AAjwIgBw4EfoPxyQ39YOPb36+5D+YCPiYvBjF4OHAj4jx5B7jyU33KP5h+QdRHjfJ+PMr02KveNrjBW6qzMfFyzN/X1U3ZXMPJtiLxHrtjHSdxa6rVVJUWz5BF0kBTkOiBgMBZXbnAxBDyIPic5AAo88FKYQKA8gH1cBwPtxzzwI9BUdHR0dAhfuJ9+PWbQPNieq1rh7SXYS00iDtOPJ28RZKnrKo5s6r5GJJfsvMXM/ZK3Ct1WCg2KShcbWl9DpqszJxb0znwJqaz3G0z2dZxM9v33F6bZKvPkB471ExLG5HrmAY+GlgKa2YwzDLJVdEdp6g7KhFIxal3pFPTQRQkzFjQCVWKkwbu86XUPafBNdyRNYEqGwWU9V7kwzVjPGc5AxcvIZJi4Yp1bLhtKWk0zErEJkYEYpewrtPjSedajDmZuzIpina8VZazxg3HlCzbjqzZB7iOm+S6syXojDH9ZhXexeMoEEinqb8lltdjhV87VxdJWYRulzuVgq9iroQ0IETAT3zt78tDJMX90/tD4botUxljHYzEFDx5RoRpW6dSahULpDVisQbAolaRsJEtKci1jmKRDeKTVskVIgF4DqQf2z3bB/vb0T/RL7/KfW7mJMz4wzjUoq84rusBeYmTiWkki4YOgTetm74qgtyy8QqBJWBcKGRWIdrIs26xjoKAUhgTEQkxN+BzG80BTRSSSOd0PHw6yqvn+4aGEQUVUS8B9Up0k+AOTwE/JvEFr/tnu2D/e3on+iX3+U+j9s92wf729E/0S+/yn0z8EyCAfSH2Dr30yf2Q6BXDjvLdsx0n6LDbGhuXhjB8O2NFX5AjlUOfFFRT8IKcEP/WASCA8e4daI9xDvHaZzmruScba67j16HzVmJs1xVim2QQX+C/Ad1t/rka3CdsZKslI1irQSbJYkhJwTOcl0weJg1hnP7zx6K3qf+7nBPkpzfSUSlKf6hAeAFM5ilOX9SCPv/n0pq9vnGxvc/xhis6wSeP9KMVq5ny9RbOoopWZbKubnbUmsV+pMUBXTOZtWMSYzyuRWbliQb6rhaUAgzyQS7/4QMXxX3Se05iXGFQxXR9ksS1WoU2ux0JGVus1y8RcHBxLRuQqyES1QpiKRClX9RQpE0UgU8/IfEfbrPle8j2xv36hNw6Eh6pC8FRgr8ksqQfLyKoYtSAfVH29NTkRJ9QAAiYempCzbnKQFEyqAQ3mQFPrBMR4HgnPPgUOPYA446/fw5REwidUeTFOACoIgQS88emA/wAH3/LoOYic7p2hGNO5hQcs0DZanHx1sdguy0zZe6v4q6FrpbVgt1Dn1lrkGs4rCa0HKyhMkZcNIsWLZ1+JQatPXOiMWl6jMEu8r2yyo8E27x8dJFQh1eYO+Jem3NyIIlIFTN6hg8R9/p5/MA+3Ugdz3GF3vWpV0tuIK3J3HPuB5WD2A16r7MyDlB5mbGpX6lS+ZQ71ZvHT0UVKWkzP4SRWTZPwKkK5uUSc7Q4dyXR84Ywxxl2gWeMulQvVYjbTXLLHFcEYzTdZumIPWRl0EXCbP4gVyAiqgkYPH2JwHQaU/tnu2D/e3on+iX3+U+j9s92wf729E/0S+/yn0z8EyCAD4h7gHXvpk/sh0Cv/ANs92wf729E/0S+/yn0ftnu2F+W29E/0S+/yn0zZ4qk0RFUScj5AUpQ48jmHkQKXyEC8iACP1GKHt9/yG3qOzolUP5N1USJqKmcKCKBEhKACCSvBDAQgF8jCryIh4iIgAfYOcPJWx3aGcXW25w1n3lgNSdgbdPvrpbLxiKCvzOi5bvzlUVY+y7KY0Y12Ajc8DFi4lflkdapNEjYZaWEjsvxqvnE4/wBJQ191hlJKibQWFrnhik4WQxflzUKtzVqsORatWD/DT9yzBjC6RmMYnCsjIKyUGtBQFXuGTY0fXmEHc6xFizNJNwve8FsyVa5rC+g1YquwGYKZMOmWVbNd7DYKLgvEbFBT4ZaOnb3GVmzP7JfpRb1XVGi63W5utzrWCshpe0wQoR5JPQW46nReXt28R4Mvl5uWzdpr1ZfZd3FytdZl7MVioVmZdwq9f1ng8NSjhxWMZ4c25coS1on6dV5Z/Ct1sDV0jyPdi0Ynbg73WjPVT2dwpi/YChsrRHUrMVIhL/V4+7QidcuEfDTzQjxmxs0Em+kUomWbJKFK4aIyD4pDicorj4gJp+6s7BIoHJ+6KX4YhkERMJTGSIAFKJW5wETC3UAheCiBPECFDw/S8dAdHR0dAdeG+w/5D/26968N9h/yH/t0CDu4r34NatBs6q6p26JsCOxU7UoWfx5O3+OCna6EXtXrmi213yzGq2WxwcEgDM5LBJQ+O7K8jVVI8pIt38Tylqow280r2VZxM1v/ANw6nWKvzZPiX2nWKG2SYXXOAhpYCrWvEuXJMKm2Da+rIroxCMHLXii0oyaLOQWCGSNKrkR3V73Hb/xxu/rMxtNjxQOTco63WiPypj+EqzgsDki2VuPFQbVjCt5CSAklj9lcuImcnJGIM7OVenxgfCOjFTOlqZV8/wC5/b/pjvZGl2WW7kvavyRM0XI1ZyW+u1qtu3WFMT29vNObRKRjSwQqg5folRIlX28S4nbvAS7ZNwqdtBFF64KmG/WMu6n2isQUap40xrsdiGjY+o0EwrNOpVRqN1hqzV4CLRBBhFQcU0p6LWNj26IFTQaNkipJkIBS8FAA6z79s92wf729E/0S+/yn1tTr5tXr/tXT3uRdfsk1/JtZZSbmvOXcH8anJRc7H+ISMRKwsu0jZOKeIHUIUCPmyBVzFV9M5wROJZ8I4dgt5LJt0mxQTKQpjfvlRU8vMVilKIIqI+JQKkmdQFPUERMAE9wW1+2e7YP97eif6Jff5T6P2z3bB/vb0T/RL7/KfTQPTJ/ZDo9Mn9kOgVw47yvbMdJ+iw2xoTl2JgFBseJvyKbg4c8pKHCoKfR4+RxASiAiQOetDe4j3iNLLPrTdcW6/biVuOy1miSgsYY8tEQS/QaFMlp50eUk7LZbCWrpyNdqraBg5aGeOINhOySyky1aliVG6zhZDoufJ/7uYEgEqgiAE8SlOAm4EQAxDGKU5B4+oBH/AB4HjpTFlRmNkO6bUKyq8TlcS6MYeVuGUMd2dyurAu9hs+uK/M63ZSoEF6TqOlbDjikUjNUC6tEoeCmK2S7CxhUZJtOSajQLXj/un9pHFdAreLsbbG4mp2PKTCsK7V6nV6zdoaBrFUi0SNo2LhI9rTU27ZsxbpIIItUU0kyplACgUCgHWZftle2KZdU//i+x+mJigUh1a5e+PFyPmYUjhVDGOfhMPWIJSeRhKbyHx6acZm3OBQUTBQCHE6fqcG9MR/qk5D6Sh+QB9v169UbJqAIHE5v3gKkETj5JHDngyQ/dMQAwgHj9gEQDoOPK69zbT/BXezpmyda3Gt1t1i2119kMQZZq9Vpdm/2U4yy/jF/UW+Hrnf1ZBWMNJ/O63JZKBjKNYJ3KwYA5bM2LxGSdLNnhId5XtmqoFcf+LSiGSK9L5qrRt/OdBNYFTELyan+ZVRAg+CJAFMoAYBVDgPKIe/xjOzWft72bL9DvN1oF71HyXQ9qqpJY/jFJKzy8ljss3Au62xO3esnkYSVibpJLLSTIXThMzNND4U6LlVVFkWs2fqTtPgPDeyOO2Vlb0bNmNKjkanN7eyZsLKjW7jFIS7JKUbsZGVatZRNBVuWSboSDpuVfgE3S5SgcQ1L/AGz3bB/vb0T/AES+/wAp9H7Z7tg/3t6J/ol9/lPpoHpk/sh0emT+yHQK/wD2z3bB/vb0T/RL7/KfQPee7YXA8bb0Tn8uYS/cc/4/+U+mbPD+gmQU0hMdVUiQGDx8UvIDD6qvkYOEy+PiIlAxvIxeCiHIhY5mwxsEyeP5VwzimDNM5nUrMO28XEMQIYpfUfyLhQqLVJUDCZJVTgn0CBxIYSgYOc7I+xvaKPdLbm3WLeuE1Gz1cLFI3e3W3EsHf22OMu5DkHKi6Ft2SxWxr0BE53fR5Hk2ESjapNEWRpmVUSdAZ2qB4UW/pOGq2r0/K4x2puzLOTxJc0bi/J2o0K+tMze65ErfCGtWWMa3ZDF7bDFltplGEhC0ysWDJMaqgEsQbGT5c3+O2a2r7veRMp5ptmgHalpb/MO4CEiMFLZ3tkOQ2p+FGkYos2vdisVzSCaWts5SZv5HBv681q7lospNqOG0isRsIn0CovaCjs6dxvBtv3hmj7W7aUaDDYvdnJikk9msH0GxrKsEsKam1zEk0m1gkaDOMZy327H1yQQaP42Pxl8KessDSYlQDrCwFmGubAYfxdnCmxllhqjl+hVvItdirdFJwlmj4m0xbaXYtLFDpvHxImaQbuyJv2Kbt2RJwVUguD+mBjTN1bWYmUUE5kiEEoCmH8HJfD6RFIxREVElBDkBMBRACl+n39rl0B0dHR0B1+TjwQw+3sUw+/29gEffj346/XX5P/Cb/wDqP/Yf19ugQR3B+/hrPovniQ1DmYexJbJSNVrc/R5vJDD8E6wuD2hmnJsi3HL0WrabTCRqLL1Enz2GxjZHLWXPHx/wJyPDOkNYG+2ujmyiEVO9wHuL1XIkHOoA8e6nY8iMi1TWyNgpcgPp/F2XoxKtqtNoYqKl04gK9Y7rU6g4TTiFHgQrdR+okhvV3fNKcbbI4vx9nmwYOqea8p6i3lhlSp1CRhY99N3Wiii+i7bj1KyvkvXqkIckpH5NfrR5X/pzWOoU6bRVZJFdC749zNnnX2mUrLddtOQO4vpnlSpxUzji549q8MrsBjeoOWjZ9RrDZH9ksUQvmylv6go+lL3lS0WKJvDSXZwqbOmzX4ikl4sL7jvuu9pPFdLquPMf7KYmp1FpNfi6tUKhV6ndYquVmuwjRJhFQ8JGNaei2jo1k0RSQas26RUUUU00yABSgHWb/tnu2D/e3on+iX3+U+t5sXZcxpmapxd7xXda5f4GWiI2UZysO75H4KVbFdsFZGOVInKQa7tvyr8vlGTN6QSnIo3IdNQpZCSfAc5/NH0kUyoh8SfxBFyqsXy8GYiPqLel4mKr5ppCAiHiBg5EoLX/AGz3bB/vb0T/AES+/wAp9H7Z7tg/3t6J/ol9/lPpoHpk/sh0emT+yHQK4cd5XtmuiAkx2xoTl0JgMi2PFX5EjgQ/iTUOFQUDwBPzOJRKIGEhQ+/A9aBdxLvD6S3DXyZw5gnbyuNMmZxs9Ux1UbRGJXuGjag2WmG1oucpb7AFZJJwdWf0iuWWpKlhIyfdv3FgZxyscmwePHbXo4fEEEBBIBBQf4eCgYBMHuAGKYxQMUfsIc+38XAiXjpUbpFfY7ukpR0gZGcxporiFm6tNJtKh14tlspnUldteH8rY5iDJumDmYqmJW+TqdIWx4eFnIhK2uoiOavmEs+XSC1Urut9p3HVJruO6Lspiur0WmQsXXq7W69XLxFQtaqsM3RZQ0VCsG9NIig1YNUGbVFqiRJMiJAKXgCgHWSrd5Dtj+K5ybiUJAV/E3poQV/SOoU3uIKCWpf+oA8D6nIj6YHTAODiPTURZtz+HqJgp6ZhMmKnBxIIjzwQRDkpQ+xQDjgPYPsHX7+HLyIidUwifzDyUEfAfcPo5/hDgRDgPbjoOZCO7rWhuMO5hLXCp7MUpHDWzGvp5fNtykom5pw0jnfDEnS6NhSs1p0vWSOYxdbFsxkJ9I19m1XRmF435ksukePIQ7KU+8r2zCJEAu3mPlEEFCisYYO+IiCBymFJEqQVM3mJQ8QE3ICPAiIB9uso7omNbLP65tsx4zgpGx5v1dyDUs94aZIqJuY5rY4ddamW2Tm4R0oiwnoiLxVc7/ILxjw/pEctm0mkVR4wbkNurj+8UfLVOpGUse2OPt9HvVaiLvSbMxFwLOwV2eYt38NMs1XCKTn4B/GvUnLdNdFFUE1U/NEhgEChoZ+2e7YP97eif6Jff5T6P2z3bB/vb0T/AES+/wAp9NA9Mn9kOj0yf2Q6BX/7Z7tg/wB7eif6Jff5T6B7z3bC4Hjbeh/b84S+8f8Av/5T+3TNHqyTRMqhi8iY4EKUADyMYfyL5CBQECgJvqEAECiAcjwA29w8Fsiuqoo2M2QRWWO7WP8ACpplJ9RirGEBKiRJIDnMt5CAATyP4hyJQ5ubzsJ2iY622bMGqG9cBp5mu1WOUu1om8TV2+p4mypkGZfLLvLlsNh5lX67CZvnQaSE2SMXs8kgrHuZRd4k4MoUSKRkj/SX9eNcX8hRNo3Y5hcoOly47yZpzATWQFr9SK84+TO7plml3qLxK1w5Zp565hnsbUK1NZHh0E3z1ka0lUatQfNnt27l6zLZJbE3b5r1PzfeabMSERmPJd/nLDR8M4jQQcqxrmMbT7OrWN7esokdGLL1avx0EpTrFCw1gWfXeLUTYoSGjb3UWFyvvdQMR2i/X3Yl3jCp/wC1bb7LOQbFI2qFXezzuFmadqcjjKacrQeM8P5plHrTYRvR61KS9egZ3CdRblYrrRjF23Bv2MNnqNlzGuPMr1ta8RtdydRqlkKAjrDUSRk+whbpAR9kimU5GklnRI+YasZJBCTZEcuStHqa6BV1gTBQx1sYLpQoiUyEiYSiICZN2QUzCA8CJBFYBEgj7l5AB8eOQD7dHQV541NVFZqokkZkqmq2OyEoCgdscokAPEQ8SCYhjkUIBRBQpuDDwHHSp9VDqaq7P5W0dlkSxNCspXWeNOXD04No01ImHCyWXcC4xq7YHEdW8f64elQDJuEV44H45JTBpCI/CqiLbulz9xTGF1lscQmweEag7suxeq9gY5YxtBwpWMVPZbhY3yWuGv7q3HcNn8LRMmkbRCtzbpqKs3w1uHF0wei3R9EPMlaKxDO9WLNmoVtg9SNh7lLyUnka+1DHMJY67mgr4ySzhPL1BWfQUNfrMzWIYtXu9gdPZimElJ8YZFQZx8Uf1hrcdR/Zk8Oba45mtc8zwce5Vay1mdRYYQy3JVw6SNsnsE3wZBJ/PQ0Cd9DqGVulboU6oE60LGwzwCvfhdpsS5QqGd8V0jLeMbMws1AyBWI+zVWwRKsgRlMsJBIDCsio+asnjZFNVM6SChm5FuQOIpEDjn8Zdwfh3PVSUoeZcZUvK9MCXQsfyK8VeFskQ2sscRQrCZasJps5bJzLH4lczGSKmRw1Oc50FymMPQTOi+bKlIYihRIcOSH55A5eAEpy8c8kOA8kMPHl78c9VCK6LhMFUFCqJm9gMUfzD7gIDwJTB+ZTAAh+YdKeWa7TaMemvBlybvProPr12FxtERNYW2Uw7GxokNVjt7ba7LFp5yrx2y8oW6W7IV1ibXEBGwZYGHnPmsiDBguFM0YmzfTo244nuUBa4mXj2U4uWMcAjLMQlgUFEbBAuiNpuCfLGbLpi2mo9i6MdsqAJiCYj0ElTzlFlEP37gip0GDVw+VBuIA5BNoiouoLYRMQPX8CGBMBOQBEePIv36WF2wmyGT8eZX3Ykkk5dxuXlifyfiqfl0yrZAZa1mOA4SxvbnRgWMwNRDvbgMbV2UnKQdfGcejFvlBkHPWa91C32hhqfK4lxxZ5ik5f2ovFT1bwzd4aRdQ4UzKGVfmZa5YpmbjVSzELBsk4R8SRkoVtIyrcq6YNo9wB1PDcnH1KruLseU/H9LrEPUa9SoOGgYWu1Fg1iq7GkZNxIMfAxjRJo0QalWBRT4cEGrcwq+RuDCPASqH2D/IOvevA5EAEfYeA5Afvz/jx7de9BaJ1u1eRbxg8BYW8igtHqFbqHScGK7RUROVBUn1JKimY4EU5L4j7+QDx0r7toP22JoHOGlUm4TgnOpmUJuBxLRpAvxdoY6ozyhi4Htc1Mk9Q9hLaTwl2KzsD9Ys1JDFOPmDJr6CHqNQdceiconFPzDx9QOPo59+fcQDj249xAPf79KUyuzca59y7CGZGTM1UxduLQXuv+b7dIiV+wlMx4/dxymotHikOVnNbk5olzzQcXMSgRlNg0L8/dt/l0b6gNoK8S/dlN5EOdPz8D8AYoAHsB+OSlEfb+sIf4+3X6K6TN6QeJgUVEQ9IePNMC8eYqcCJS+ACUTB5CPuHHPv1aROmqdYDqeKSTcqrsBTL6DgFQN4CCpxKYgIgkoKgGKUvBwE35cLevO6VuypZZ/DWglfpmbstViQdNsqW+7yVro+DsXsWxwQIEjeYysTD+x2yyHF4tRW1Lg7NXJX8PzhLJYK+QsYMiG2+xGzuCtaKo2s+cbvG1CHlZJKGhm7hpIy7+wyyxFToRsbDwrKSkXBlDplSUcnapsGqyrYjx239dITaPx2N9p95FHz7O7q26i65FefijG+OsS5Rstf2JyE3lhMeJPnibhWsOji9SsRya7GbxZULNkap25WxKGk51v8AhuP+O2C130wreGrhL5pvV0u+dNhLHFqJzeS8gSz56hAuZw5HVwZYgobuQkq3hSrWt81jHM3UKCuxhJIYSD+MRV+TsfT3AFoo6OQrghQIiCKzNymoJXIqFA/xLZ6kX938PyKQlRKqqmqIcqlKKZAEIItM3hHTzXmWm0mFRxFhPA9DOLGLYMgr9PqMRDIEaxUWixiGRyRrRNyduxZkj2SxURcGEhCgA9QN27MR3qGxTN5+zjU3VU2L2pmU8s5TrlkK2lbbi9hJA4eVHALq1+qu8s9Xw0nKTMbUnToWqaaE1Ii3jWALqEUirZ5452r2uxLpdFsBmsQ41+WbA7ZukhLLxqKEaYjfC+Bsh1WWFtEztQz40kL5Ou1ETTDdivi5r8ay812pumyNSJpIETSKBUiABEwD7eBQAC8ce3AAHAce3AB0Hxbsyt1RFMiRCcCP0JlKJjHH354/hAvHsAcgbkftwHNd0dHQHR0dHQHXggAgID9hAQH/ACHr3o6CzrRSLls4YuEEFY9yiuzWYHKU7dZosmZISGKYviUFEjqJrJgQSqFOIGEQDgVX6nvFNWdiMu6DSvEVSJor3O+naLlyeNYuMaSzoxMmYXxrXGJHUdXcb6wrusfwbIya8cV0nfmoMoNumgsUG0dLi7ieMLg+oNd2QwrUXNi2G1WsbDJ1EiIYrGJnstVJsVULzgR/bxcN5GIx/fCkiJu3synWZSDqmwQuY56o2bmQDV7a/tZTJLxkrbXt25WndYN1rhZaderDHtbLNw+t+wtjpIT5fk+ecaQpTwU2exnsjtWWuj+GsUwgqnySPefEqmTuGtfdUYqZgtuqPcBolf1I20g8hRVap1fLKzFlw1miKsScyep3PEmVpOAgPmSUqhCSS76Ps0JVHkYVxHkbIu/iXANWkYkyhSc443omX8X2mPtmNb/W4u206zxK0imwn4OXblcNpBAHrVo6IgcviDdJdBMximOKiaYgADGW1+qeFd18NzmAtiqQlbMfWJw1lTpFV+Xy1ZscOc4w9qqFibcSlatTIrp4nFT8WdrINW7h6mm5RK4OVQNqiOk1DHKT6hIbx9hKIG+/uXgRHxDj7iAAPPsI+/X1SVIqUTEHngxiGD8ynKPBij+QiA/fjkP0Eeucifm99e0AF8tp21h3v7fCmUq3YnstZLta7luhrtj+3J2Fe2lRbTjV8hluj0x6nXo6OVmMhEnG7U4KMYEAcuypO11u2m182xpT++675SgMoVWIsEhV5d9EN5SLeQ1gjFTIvYmYgLBHQ89EO0lE1SEJIRbYHApLC3MqCSglCZLbLM4CtTc/IJuFWMDFyE07Iz8fjRbRTNd84+CAyiRTOhRQOVIplUimERAVCgIj0tPtbFJkvDV83NkQLLK7wZTn8+Y1nJUnrXxLW+zO3k1rpQbouf1vhH+OKfPyESyrbWSlYSrg9esoaRcNlzqHyvun2qyI6xNML4/sUtT8rbY5TpOs2J7hGP3MS2q94uaU3aySlgmI9QJaJrylcotii5B3ENZF8Y8i2agwVbOXKiO7tRqVTxvTq7RMf1OGo9PpMZE12rVKpRrOFrdeiI5uLVhBV6GYptY9hGR7dMqCDFJBq0SSIQqYABSlAJG6Ojo6DF7rDEslTsVbUcKNCWOFlYE7pFMFlmxJePcsTuEUTGTIqoiVcVSpHVSIp4CUyheeekE9guwExdjLant0yl3ylkyxdujZSzYZa3XIsSWIaWLFEnJ2FriM1JSNNy6rWBbQVKkAJX0SoQUIRdq2h3j9A5lE+hV2YpUuDFOYpzeBhJx9AGKb6zgIgIphxwYAAwjz/CIc9c1mfZh5pz3/ADWDI7m4ZWmcddxfBU3r7LUGrR6iGMavlrE7ukkxlY7WJJRrFr/GwElfPh5By0WnIgBXawbJ40kJRVuHSj8YkAlIYDEUMQTimf2MXj+qbjkoGHkPbyEPv7/r4L1IEjKGKYolABMkPj6hRH7kMHl4lMX8wEwByHACPWOS01Dw0fIzFgmI6DhIxom7l5KXcNI+HSTW903TmTfKoINEERKYpjLmST8lCAJufEB58Nnu7hlDNea7hoB2pqA5yjtOs/8AlM7svaq20cakYkaQyi7a/wBkmrZ6cwva56hzikFBr1uPqkqxkhnlnqDtYkeU5ga3u3v5rZozh+fynnOyLrBFu2MdBY3qLZhZ8p3qwSaTtzDwNNpPx7VzIykqyYSL1iL5eNYnbsXAqPUlPSIqjiOxxv8A98h+8n88P3WqXZxycsrP1vBqC76s7YZ4qEWf06j/ALRJKJbpfg6g5QipZSyytdQtVg+SyUHENCxr0vLltufpD2dqhjbMsBvruxY1Nm+5LLxr+Rt2VXy71zjGlSM0qyejAYZpUmBIuvw1DWScRFHsLWHr8wwg371q2jo9J67b9OzIcwuSpEAirYEVTKrFKCbhdyB0/TK08Q4ImgUVSqCY6fgJiAmBwMYShpTj3GGrfa303XrGMoiIxJgnAtVdyy6z9ZYvxswRFuwJJXCVatFXs1ZbjMKR6M7PLoO3slIqkdOzqmKJwoe3jinIsJiGXzlnmsP6lsbtXOhmvMFRsoN5OzYlWsBF5Gt67KWQF13FiqmBkJWUqNLcLC0RRjXDkGsZHEWOgaHtlXchtft1ijTiMijTGEsSFjtg9r5RISTEWWTjyJtMI69ZEqkuZtGTdVzfFTt4u6zpAZlswfYiYA8aFXXamK2dAABMALyJeR8TCPkJwD28xEfcRN9xEfcfz9+gp27MrdQRTKkQnA8eCYFMYTiAiIiH8IF44AA5A3PI8CHVd0dHQHR0dHQHXggBgEohyBgEBD9QEOBD/l170dBYnsI1kmD6JkGjV3EyLN3GvYtdMqrNzHukTtlW6qRyiQUnDZRRBwh4mIomoYphEPYVbabnDWPPWX9C7A7M2iI5V5mrUX5oYUk5TA9kfgvfMWY8rqAuY6rYx1cl5qgY0hkm7toDxpNRZmkGzRROki2bpa/cTxrak6tTdqcOVlzMZ31Ss7C+w8PDEZREvl7GbhNeOyHhWdtpV0JNhjmQRex+RbFGB8WykJnHNdMpGOnCDZZuH7yBos3rV0seZNKrfA6m5uuUvNzuT3tZx/DzdBzu5kngyjpXKlDM7hoWcvSkkQjeKyxIhK2ipxklZGsS2cJWB+Q14wnuUzstkHE20OPJvWjOtejTuiR10dxhcWZNdwThrD3C14CvISHr2unRM3Ix7KPf2+EoltkmU2wco1UCnkSMdrsc5Cp2aMd03KmLbWxtmN8g1iGuNFtMUo+SZWasWBohJxUw3O+bM3rRB2yWQcNQUbpr+ioJFUkhExerBmXAmGNgqsSk5kxlT8nVdlMfiWMa2ytxEynB21BB2zaWWBLKN1girNHJyT48bOtPQfx6qp1Wjkhx5EJzK9QN9zeID5Dyb29g/PjjkAEPcvkACIfl19kV0nCZVUFCqJm+xijyH+ID+YCH5lEAEB9hAB6U4o+2p0W9EromTN6dc1U1YOJjoWIrZtkcJRkGciVaUsEzZLFGp5vqDatpya97ybcLiGRzTLCFNGVidGdlnDBieHsuYqzNToy5YkulculdmI+NnCvIJyAOUkZ5sL9irMw7lNrNQT98j5rGYT0dHShDEVI5apqpKkKGWXacY1eqWG0SablWOrMLK2F8Rj4i/FnCsHEi6BiU6iJDujINlCpEMqkUwmEplCAIj0uXtdMD3TBlm2zm/SmJbdPJlq2DptikCi4u6+BLnLyU7rlTbo8UKdRKUxzjGbZ1lvAoP5OIq6Sa8RCyLpiUip7v3TrbPhrrX8D0KwStTytt3lql64Yrtsa+cRbGuWubTlr8/dWSWj1fm0ZX39Nx7a4N6vFM5J2urKt2CjE7J27VR3qq1XrGP6pW6TQ6vFU2pVCNiq5V6rVY9rEV2Dh4tmVmxg6/EMk2rBnFxjVIjdBimi2aoIpEBEgAmQOgkLo6OjoMTvEFE2uqT9SnCOFIm2Q0tWpBNm4UZvFWM1HuWD1Nq7S+tqsZouuBHBRAyXucoiYAAVvdr6flKjiXImq9yOlCWnUHJ1oxNVqO5bJhYoDXBrNyjPVuVnZBuZUk25tmI4eJl0LC4XNLTpBUkphszeLKIA0R2IeiYgmEvn9PkHHIcfVx7iAfWAeHuIB9XSmbM5ca3d0OqWM6AU/EO9GKj0y42GTKV8hbNrsUEgWWC6PCCT4h5XHp8KMMtya7VkihCzBIZZ7KPCySDFFYGyA8S5AogYpxJ5iQ38RR45AhhDkoHHkA48uOR+/t16V0mYUigU3moAiJPbySAvscVBARKHgPBTABhEDCAcCHIhZ1l25E3i7hwm3aNG3rvDrlTI0N+6+IF0dyoJQSRbJpnOsJwKQpQExhDw6Wvaty8hZznJnF3b0r9Ky3aKvJvI/MOVchvbZRsOY09B0aNNHws6zqsvJ3jIkgqqawVCOhYJ/QJ+EhppSZukWsMa2kQ202M2qwRrNXo2YzTcmteLYZEIeqwRGMhNTlvnTpLrM4aGiIho+cnXdrtvhUXb4rKKRfKNkXki1FZMw6WRWHtoN2xdSmzzux6uYMZvj2DHGG8M5YtMZljIEdOnO5Yq7HTsfH19OrLRVfcu61a8JQEpkGiTT6YXdL2hb8PR5n2x2uemlRwNY7DlOet18zTnS2xwt7JlPJU7IzkgxNMLISttiMZViUkJKIw7R7LYG7Wae47oztrWCOY6JAED/KWRibZmaKOlAFwTxBuCSzNwiqJVlABIxXKD9EOExSKcwemiU6yRjFIofxOmQOggPJlxxNqLrpYLU0i6vjzGOE6G2i6pXGSH4dqkMdmkzrNHp8e2imapYhB1LuYarRBGDJVFmV6kBSpopjxD3b0wdc8fYYdZPzJFPmGwmysuOaMzNLOg2e3GkyNrFzN1/BsxZCuXTqzw2Ao6ac43qMk6WSAYWPAzePjU1RalhfYJNztzuXjHVpoirIYR14CMz3syogBZmCmroowKzw3r1kery5msbJQN7gbPKZfaSKBJpuwmsXRIrIJPDN1Um0IcAmAF58f6oiPImD8jciI8+X35/P79BbhhmgiI+ml7iI+yBOPf8A9+jq7dHQHVnkEU1TiVUqaiSqZklkVjAYirdQgkXSBE3JFCqEN4KgcAIICACI89XjqkXZN3J01FiCYyZgMUfIQ9g+5B4/iTU9vUTHkiniXzAfEvAKk1bfm1V2by3pFPppw2Mr8tLbDafprCCTV3Unzj0sv4MoVTjQdQlSoOvHlj80SmZaITkhyMt8piz/AATsU2lFKYDnBcpzcIl9UxwIKLog88giiBhHlL2A/JClDzLwI8iAaGdxbE1zm8UQeesHVRey7F6wWZhlvGUDBnj4WeyYwggUWsmDX9xWcx8hBUHJBE41W4Nmr74aR/DkT8WxeC1QBLZ/DeWqLnnFePsyY3sTS00m+VJpZqzZ2RJBmwloqUSTMLtmk7atHqLYyyKiRfiWrdwQEuSo8GEeglISkbqnOBEuBQKR0YDGIcoI+QN/hkClEolD1FQOPJRD6OAHkRDQ3J+jkYpcJTNOq9yitS9hrDKSTjIOSaVQ4Kdh8uoyworSaOWaE+NEwF+sLdZskFTuFhCUm6WSRsQ186Qz0kVZg5GqIAkPiJjJp+JTmETn4EA+6huTiIfkIjz+f5j1jk+YscydvP3qjePZLvHICYFFfTaIKr+o2Oc/km9ApFARMIpkUE3CipQKA9Bz4Re14Zq7luC8PbX/ACPBVp0+xFerzlShXuV/+SuVsnZAeVZDCmcsFWNZNSEkC0saXe/w/JZPbY1yNFBZ1i1OuPSOpsWfQmRwZVob01BMJjpiU5jHRTBJYPIhEzNgW+sAKAlUSEScGDxP79J/0Mw5hjbvXjJuyedcbUXMFT3PzRYdg6Cyy1Voi+Wep4adiUmGaJZi2NjLt2ctjIzq2ngIuAfSsFWT2B8euvinfPRNnykftlog3K4pSOVN6NcyLDU6ziVsjTCZ9wrGMR8qy/j75cbNDK5mqx2y0iW72fJN1/HDAI2vkr8TMfHygNAbSX+Ev+Qf9uveoQxBnzGGdqpH23FF4grixexLGTctEFvhJyA+OBX4Ztaay5Tb2OsP1joOUwYTsVHugO1V4R4KI9SyaQAiYrLnI1SQTM4d+t4iokQociQQIJycAAGFQxTCYgeI8cD0H2k//shx9H4jgQEUBMUoKgAD9JgOYpDlH8yGHg3HuHt0nrvTZRqWKNLZvKDqSVJm7DtphM5arQkXGStonrBnzHBHqtK+HokSxkVbJFMyS8iE21mGqdUQM5jwl36HqthPMF83YuWXrlLYY0Gr9HzLkiqyD5PJd5yee41TBGMmbE5EWjVe5xFblHtzsdtMo+Upv4Di7bWDhXJgLJNwYLRXzHOcC6W0nEFknct5Gtl5znnvIcY5jrVlG+2OZmo2DfzxDhao/DFKkZB9W8DVu1mBoE1WsZta9AypYSFGSbnNFMARDULDT3NPdXxBRMy5Etk1gPTu/wBbiMm4nreBsm3ui5+ybGzyHMJI5stNdZ151jNWtA1fJyWPaHaLvWrf8+VJOOiBBRvrt+xpRqljKj13H1Br8PVKTUI9CDq1br8UxhIWGhmZClasY2Ljk0mjRsnyYwETSSETnOYxAEwiK3+24shieIz9pHIuDVtzqfliXi8M0F6J3s5AagWJRUmvVjlJ1L4oLJ+KTwV8BGTkpJ5alvlCnz9u1D4D1mptwACH4AgAKhhDwAQAfYPqEBAPqH8+OQ/QR6DxybxIQfExg9QvIlEoAQODfWoJjFD0w9gN9xERDgo/lBuacv0vAuL79lvIk/HVSp4/rMhaZ+YepyirSNQQAiTAiiEcweKrpPXi6TVcrdFZUDKpHFIxSGEs3vimMgBS+AmFQvBFQMZJTgDD6aoFKYfAeORECj7gHsPSo9lpB5tdt3h3TOJZHlsR4oGK2O2vl2KoOI5BRoYrbCOA73U5U7aOs9FzqzfZBm5MRbS6cc5xgwCQZoqOm3mEl9vrFWRK5iuw52zjWH9V2J2isauYMm1mUcR87P4pbTYqPavgJCzpvHSs7ScOpyEtGU9H4orFo3mH/wAGxZguoQ7FWoFBIAKUxAARDwN4h4/b2ApBMUoB+RQH2+3X5RaN0U00yJgBEyFIQB9xKUpQKAAIhz7AAB/7B7ewdfdNMiRAIQOAD39x5ERH7mMI+4mH8xH3HoP30dHR0B0dHR0B0dHR0B1ZJNsg5Moi5TSXbuEFGzlu4EDpLtFkxTcoCgblNVNYhvBYqgAQSiBTCIG46vfVIuxbuVE1FiCYyZgMUfIQDgPumYA/iSP7ComP0KCUgmARKXgFTapSKOq+yGXdGLAdOHoVrcTGwunjJUBSQf43kXpSZXwzRKrGldQVPx9rY6kccw0G1FxEFkEb0AxMSKbJ36LSClN5nKuBzD6RfUMp4GQdgbkBFFEDHEPTHjzDwKUvmUAEeR40H7i2KLdI4xrmw+EqorYdhtXbTH5UoENCDHws7kuvRZFy27CMrcVHDCShceXhI8dNWxi1embSbqnwYuY92o1bCjtTiHK9AztjCgZlxhZmluoWQajH2qn21knIMmM9X5xJu5ayTJF21aPW7V0UpTJldNGzpMAABQJ5GAAk0SggqoPppeJkiEdGExyqFOiAkbi3QIQSin4nV8xAQEPpAAEBHhKm2Xa8kU73e9t+37k+waxboTVzp2QbURnbLTF62ZvsNSRnmZYfOWKIYi1anTWNOyyTyRualVtFgbvEfIjVcHrhUjwCtECiQfETGTJ6ZTmETKePAByJx5MJvYPqEeesRt0klWoeSsSpXC7GsQkpLu0SCmsqohGsVX6p0AWUADyPotVE2wrGSTMKhxOsT7iHMpU94rBknul4lxv3DSY41kuetTbJdA/2KObTe7rhjLeYMuTlUldX854fs7unpVJrOxFFoGXGaE5k0cc2+EbXI7GvNHzeVnvg+n1qd2LkEl0+UzE81SnIl++OYQEiKSwCYy3wYAYigremJxUKYoH4MIJF1M1dwT3ENHsj5F2LprbIVA7gWVrFtHTncw5eN8s0jEl4kV7RgCqu7q1A83VLxh+szsrA19SoTrqMpqEtNR1TmgYSLkV4nsxt7+zk1yBcItSy74dvUMsQNlVh7JdLnkHdPXijWj8QK3U0W+tZXY5bo1df/hxnHuJ7ID22tGxUzRsCBXskCYdG/R1rRrZttgTbmpzF619yRF5DrMFY39Um3DNlLQ0nAWOKVVRkYKdrtljYayQ79qsmokYJSIaormRVFms4KkqJJpn7VG1iGlpyxSsXX4iIbg8fTMxIR0XGs2YGApnLp9IuG7FsmmYxSrqul0kEzKEKCgmMUBDIJHw+GEFDCmQxigdUBHySLwIioUC/UJgAOOC+/BhH7APXO7/SJLxjGL0eRftNiXuJs6YjznirNeFa3i90vIZfvV4rjC3IwdGiIaBcJ2yEaWWrS9offPiNkY5mSIKzO5Ik/AR/OzXdyy3nvN9m7f8A2p6LIZH2PGSUj5Xbyz16Pf6iYWZRB12l7lnFkVLKq3axUWaXg4JaDj6tNsJQs04fMVnyLIViSBqB2V6Jj+6Ptx9xLQrsv3Hb7VbVHXzNsi9kXmOqK/uyCIu4TC2PpL0YGswFCWKtHUGUja7X5ONhnbli3axrd87QMGn1bx3vJ304+NvWc7Enrj2d8uJo3Oo4FgpMIjZrP1JaE4pElkKzU/1iVTHWSoeYWsUtW427zQRUnFQrcIZYqYLt+ijVHV3CGm2CqPr5r1SmVFxjRY9NpFRiHgvIyLoUUU3lgsstx8VYrRMGQTcTlhklXEnLOwFy9cKqmE3SYOwM6Tw7jrbrtzSU1mC7T/b22bsGLozIuTmaTA9ixFb5GzhiclVRPLPXDKuxcHRHqRoqPRQr0eV20RrhnrMVlEeiFpx6JQKkRACiJQTTAAIXx4D6QL7AX9A4AQ/MAHoB0bwIQwlMYAUDkSiUCkDxP9SgmMX93zwU3HkbkxeCj9wgnN+YqfgPFWQMuX2ajq1WqFWntjlpJ6SUUaM1BOi0hmRko9g7UURk5N40j1/QSVORRdFU6YkTOoScn4GMgAF9MTeoXgipRMkpwBh8FfEpxAg8ciIFEQEodKj2NfSm2O42IdPoqNUlsJ4XLG7HbV2BksVzGFn0UyM8Ga7X2pSqjWPstIzVDTt+urtf4WYTjJHEsSD5s3cOWvISj2+sU5IqeIZrNGdqvI1HY7aCxK5py5U5VzHzs/ihzZDOZSA19LZ0nrs9gpuDU5eWqtLD4r4FrHuXIMWbMixkzMQagBUQApTEABEPE3j9PHHsBSGMUpQ+wFAfp+3t16RsiQoF8fLgOPI4+R/twIiYeTCYfuYwiIiPuPX0TTIkQCEDgA9/ceRER+5jCPuJh/MR9x6D99HR0dAdHR0dAdHR0dAdWCYYtJBNyyfoNnrB+0XYyDF4Uizd5HukDoPGijVUDIrt3KSgpOE1eEzlN4H5Awh1f+qRdi3cKJqqkEx0jAJR8hABDgQFMwB/EkYeDHSNyQxiEMICJCiAKq0+fp6w53zFoZZn3wdcVeTWwGorR55AMzhKwS6at/xrTa2x+KgqdjzWOcsdAx1WIwruNF3GzLQ8RDkZtHAN2iFKbzUKuBzcJkBUVPAyTsB9hOgj5HECpm4AwCQoAJgABN0vvuL4ttI0al7Q4Zri8xnfVi2sMjQMRCGYwsvk6gERdxd7xDY7eZwxlI/HEgykm+QJ6HbOnDWUmcfV0Vox04btjobg4xyXQs047o2XsY2Vtbsf5Bp8Rb6Vb2aT9iysdUszVpJxM0ybumzR41aSbJVu6RTdNWztIpyEO3TETlKEgCUqCipgIl4GRIRycROByqpFBNv6CBSCUUxKc4HEBAQDxDxEBHjQLI+jyMVbJbM2n92h9Tc42KRmH2S52p0iGmKTm5STcjLyo5SoLoI+BsN0XmGjZGGyhIM5i31GHf2ZlBHFCwyKK7FAaogJB4ETET9MFDGEynjwAcicfqE3sH1CPPPv+Y9YRd5dvUq3YLO6I5cRtXgJSwP0khTVeGbQUa5klfl5llSeMgog1VSRBRRBFUyhhVXLwAiCCaLtKhn3ufY/p+0RIXBFl06wFZhvuHsiy4qYzt+f8kzlCc40z9r3YnCa9Zl4+Iq8XkGEqk9fFMeZiSgr2o2jKSVlIWkkZ0IeqoZv4cnKoC5SH9X1EAIUwicEk/hyqeKifh4iYo+mIAYAOJRERT1oVgTEG0+p89nbYDH9Jy/D7mZXtW11Ic5Ir0fcrlUsX5Feys/r/Wpl/PsZBeKumHqNbXNbggrj+QjKYZZ42qMv8CPrHz4z/a/QxBui3ZZY3o1xUUCtV2LYN6eGedfoSBASQz2w2Sz2SIc5porWrM5BxcMg2y0zWXJGwN4QGddmCy0u5bA2zo6iDF2bMb5pqzK4YnvNevsG+i4uT+Ii3ZSvowks1I8YNbBCKkQnazLOWxjqmhrFGRko1VTVbumbdVFVMkjryiTVBV68cIsWTJud2/VcnTJ8OkQnmf1j8mSTRQLyZwt6nil48nMBPIwB9pT2bc+l6/icB9ITFADcfmYpxAigF/jAhvuJQEPcA6Sx3s8zY+wVrDWM5SUg8Vz9r7kys5v1jrdfipe3WGcv0T8dULC6CjRjJ43lodviq35BRdOLQRjVY2XWjFF5NKQJGged7ZudkfYC0yuLO37BUfJMxU30kzyvmHLCdwquH8cC0dDFJRFblGtYknuQ73KgseyVE0DDTOOJGIg5UJi4Ryy8W2kZCwto5izEid7nbrYMgZryxmiuy9ayDlvJFrnrHLGStrY57zAYtj5WTk2WDaRapEwzLjHmNzwVLauY+IBBl5Q0X6QayY6rOX+5nj6mZlzda5LDGqVzgoHKeEsZa+ZSvtXyBkmpXxkhPVCd2HtUQxrDyHcM6jIng7XhevSd5x7PO5ty5kJZyaCjlF230WrV2j1Gv0ynwkTWqlVItlXqvXICKYwcFA1+IbpsYmIhoiNIkyjoyPZIotWLNskik3bpppESTKUCgtntjSs1UcVZP1OuqwwFs0+y1aMXVCirpkUlqlrCnOS7XUtzLyjP4hvYHFpwxEREoEyu/f2F94nc2r4WXVVSM0tuAAQwFKQoeopwBAEA9zc8jyAfUPPJuOQ5+wiHHQfh0bxKmPgYweYgIgJQIT6DcGVExi/uxH6BAPIfIxfpEOeNf895mq2v2Ich5buEnHw8VRa6tLnUehJ/CKzUg5aw9TgjlYsXZzJWGwykVCGVSTU9Bw+ScrFIikosnPj8DCkXx8BEDgPpqgYUlAKAmEqviU5gKAAJwECj9RSh+YdKkzyDzbrdjFmr7Mjl9hPWtON2C2SkWKhXcHN36QjCs8Na75DqUoo3YTtTvFbtc3lpB+VlMJxk9jCC+ISavjIGIEt9v/C17xphl7kbMsdLsdg9jLA7zXl9jOLMZizUOVubh1Pw+DHVkI+eKWGrYGZTbrG9KOd4do3gY5IrFqybmBAGDNuARKBSmKACIeJvHkvHtxwUxilAPsBQH6Q9vb7delbIlKBQLzwUC+Rh8jcfqJh+oRH7iIjyI+4+/X0TTKmUCEDgA9/1ERH7iI/cRH8xH3EfcffoP30dHR0B0dHR0FBJETUbCRXxMiY3isgdEi5HSJynIo2ORQDF8FiGEphEP4eQ+xhDpVOrzgmru1eWtJJonyrGmSlZjZLURVdVVvHvYR2uVLNOEaFVYb4qGqOONeAUx2aFYvE6+3kxyM5+TxzkGT30GtvDeLc4gQ6oh9k0xAFFPv8ASQwiUCGH8jCYvH9oBHpd/cCxXdrLiSvZ2wfVHln2G1dtUdmHEMBFuGFelMlEggUWsGFpq3uHkXJRmOMkETj1LtG/MkmksNdiPjmbsWjf0gYaD1ECgJuSh5gQohwcBJ/VW/difxTNwPj5eI+w+3S3u6bbbS01jcYlx3NzVKyttjfanqfiXIMRKOohvQ8hZa+YhX7jYZCKcEm2lfjU4J6V04gm0hMplcAVqwVBRTjbzC2X6NnbENHzVjWdbWSj3+uR9or8vFsHTZtLQ7pIBKmi2k2bF8iQivrJlFw2QV4KIgAgYR60FtC58/8AdYqNTbrAFP0PwgrfMj1O0Ao7rluuezjtAMJ3WoQxCvo1a34lLhq+lGxzjaMlq8FxRCovXASM16IMdx/j6p4yoVTxxVYKDrNbptdjYSOi6zGNIODiGzJApCIRLBig1bR7b1vXURRaopFSE5hEhDHHnLRH6jesJWhilVAiBEQVKoiHh4uTcAYFDAIjwmPI+4j4+/VM5EFTCVQVWaJmwLuDrKJqJMwR8jFM6IqoZBQfExzGKn6qRwKIuRDxS5Xbbty7tmG4zGHdBYOl5evdPkZL/adlHJSeQaxgbHMaxMijHw6FuZ1/4/IFsuZjyS1PeYxZXajphWpUtrsMIV1CBJBGG8WPtf8AXedR2NoeYUNLNgLlOvm8jbMW4n/2jymf15b0FpJrkrD9bg5cuRpBFy3ZpR+SLVAyzzHoS8h8FLxQWR38Zo7ifOO126s6fFndox7Y+2Zjps2CxYqjccZufUdfZmTYKAnJSs5nOhWSPfYxSpCK8UYcbv7JCt8nFuqyakRMhUTAwcDrrpPR8IWyw5osdluucdiL6xI5veVsjWCZniNZZ0ZVWaLiGjTUg7p2Ca7YljohOVPE8bU6/KkiYb5jHLDGMAb7JZLxBjTM9ZUpeX8fUvJVXeOUZD5BeqrCXWvNJtqRUjKVRibMwkWacoxBwsMXIkbFcxp1FjtV25lDCYLvQqVSsdU+FpmPoWFqVLrkJHxlcg61GtIuCioxumcrdGNg45FBizQ8R8/Fu2TKuJhOIGH36yRBBmVdQxHJlDOGqJXDYTGBEUkvP0lCsufBmZXyU8zlTIZfwADCb0gAqs3dd2x0DZlc44Z5S3r13TULXYbB7FSiJZ9wxCNjHPXXkJlDI9grDvLtdbt1JEt2lsn5BkrozBlXy1ZlJi9lgQ3nwznHFWwNZbWfEN3gbAgtHxb+eYprAhc4VN6VwLOHuFZkE21pqbsTovUwZWKMjJBM6C3w6AlBQeg0VyumOu/ciwpnBgRSs4t28prrX7ONulDBNtJnL1QdR6mo9Nh41wZ3MVV3MEtWaB+IgI1rDy3wxQuL5qLCDBRrxHiZE1RIBlvRUEp/Evhx+olA/gUShx/VEQ9+A5HpfHczxZf8jaoXaw4WrDmz7C4OlIfYHXBFJ1HIC0zdjoHqtOerNZd40g5pvHJScoc8HaDKV596vD5BUyaQBtNg7K9Lz/h/HGX8c2NK30LJVTiLZXrGzaOo1jIMHjcOXKbZ+1j5JBJ05KqAILM0RS9EeEiAJfIKnO+aqRgjD2QcxZAnW9ap2Pa48sc/MP0JJVszYolBJuqulFM3soZBw+WaNTiwZOHKQrlMZICeRg1V7emLLzU8VS+ac3wCtd2N2gsjnMmXq/KHjpmexqpZPJ9CYPZ29u6fu7BQ8Qg8lI6lInfqNo9tMyAM2TMHKpVIk2QXebU7c4j03YNUZXF2HXMZsjtk/j13Tpu3dAIIYNwXdK7Pg1irXjrOjJzkSbmmiTabaRrrGUX8xbM1HDMFmpJtjtvR9IDAmiX0USC3bmU/qgAEVEBOUBKXgfI5Q+koCPIAHQXsPsH+Qde9eB9g5+/Ac9e9AdHR0dAdHR0dAdHR0dAdHR0dBbpRJFdodFwUirZXySctlG6blN43VTUSWaKpKlMQU1yHMVTkBASeRRAQMICrHVJVDWTZ/MWjkyPyuhXBSa2Z1FIsoq3YP6LIySDbMWIKLWIcHUFT8aa3P5nGkLXYt0WASfo3zyhYldGPfGatUej4tzm8DqiHuCaXAKH+/wBJDCJfE3HI+XkUfYQAffgV09wPGlvlMY0/YzC9Xd2HPGq1xjcr49hoxdlX5HI9aaN3jG64hsdsVdxkq0xpamLxnarVBi/BpMStErJ3cc8XYs/RBioPUePIREAFTxKIcHASD/Atyn5eKZw9y+XA8fcA9+ls90i0WUNeIrC1BmJinZO26yvS9WsZXqOk3MYwpttuqM3biWKzu4twEyjXz16gT0SupDNJGSMaVSbgyM1WdnJuZiPLWPc5YppGasUWVvb8a5GrURc6XYoyOdsm9kqUy0TdwzxFnKMo6QaJuGjhJYpHTNq5IUfEyZRAxQXtJR77YDuxNVkH6ZsfaD4JUa3yh2gzh3DWjKm0Lqu2jDWT6LBFTfQxrRiutYoyPXnFrmEoqxwJL8dlVHLmOmZ4xAZhUqXTcbVSu0Oi1yFqNQp0PHwNarlYjWVfgISJj24NoyFh4iNSasGMe0bEFNpHMUE2rZJICJpJlKUOr0ZukgYh1noGTbpKlUIcgKLLFMYgeRhTBRdYW4lAhicGKf1OVfcpebNYbDDVmMe2KyS0fVa9EMU5OXnbFIRUfCwaaQlTIvJSky5QjY/2WOUXQOCI+5gWXA5kgNz9bHd2rLuxufLdoN2mccyN3zy2l1Yub3fslfaSmoeFmcKsuzyC/NPukZFS92mnzisDDOYODq9phZcJReQizybePB0iETd5DHeu+lUZkje3WLOlk1+7gLq544tUHhjE9msEtDbW3qPYWoYym5H1drDlxXrMveo59PTL/IVlpj10m4ixO4kwcyBAWWlB7YbJ92rNtXwn3aLSn25NJ8ptMOzFV1WiIHJEHYNm5K616xyVchVtiYWpM29To1gYN3sndqRfL9WJptNNqqklVHijR2eNf7o92dqHhLLTXeTauYdbOdx+ytbDIZDzlMys1K0OqyNleR0k6quH8e2Azas06CqLpqtH0OZhKtATMTBLPmDUI9F84bKss2a1cwjuLhqewfnmisrxQLWLU7hl8U/iJaEeNjGVY2WoWJgZnO1CzxoiolE2GAeRk3GoO3qDV23SduCHD86taz4Q0/whSsBa80qHoOK6FEN2EFDskSKPXKrZFFutNz0s4IWQsVjlATItNWOxLKTsm5D4h+usudQ3U7pCmeRKoU3gLdsdAASBUGzQHAkOUgkAgEdKuARExVCgoVEEzFKcPU4Nz2WWG3y7OkHfLnTkJ7fHt8Dliv2d9RrVa7rkPdXX+h2kLAvfHMFOWs7s2Y6jWZIK0yZNbRfrBbCNjEWjIQ/qyIkcrrZtnr5t7U526a9ZNjr5XKvZ3dTsJW0ROVecr0/GKOEXlan6hboeu3GvSCCyCqaYy8Ex+KFu4+XKuE0XJiAkzYqZk9Pf6QBqxl48zmix4/7i+GrFrrZ6lXkTJYoo2UMaO6WOK7TbyqyjSIkXEhXnd9CPcSLVa0RKYvEa4yWYvZs6XSQ2cFRa+XpHMJBAT+JxMJvIB5VAywlOJR8ffzAFPf6igPt0k/v7Yfmsl9uy+5Dp9myHVr1qzeqHtPU2+LkXg3qbm8bfOYY9YbPohRGYZpSMXcZRU7qKcg6IVmCZTA0VdCDIdT9ga1tnrRg/ZSpRlggqjnjGlSyvAw9r+WoT8LF2+IQl2MbYAiHz9gWT9J0UHaUY9ex6aqJytnByeIiGaZ/zZTMDYayFmG8zTav1WgV9zOysk/SkTtkiCJGkWVdOKaPZT0H8u6jmBjsma66BnZDqkImBzl1h7eeKb9SsQyOW86QC9Z2R2csTrNebKzKKRs1PYzl7adeXjsGIW9q6kHNko+EjScpWKL5yC7djGvHJWTVmm4UIeHNh1pPa7cTE+obFghK4fwYtH7HbYTEc5cu26VnOkDbBOvt2rc/8NG2jHmaoWYyFdJNFFpNs4qUxPBA/SYuF2YKtWKgqkZM5RMJUwMmn5IIHWMYwlETAsYBU+vxETiY4eQgAmER46C9dHR0dAdHR0dAdHR0dAdHR0dAdHR0dBaZpq1esFmb9JF0wdprNH7By1SdtpFk6QVbumLhFYp0zIOElDFWKcpinS8kzAJTmAVa6frNdaNh8yaHWJ8q2rbx3O7IakFfKuU0Z3EVhmElsm46p1fivi4Gn401lsFooeOqjBOTwfxMXNs1IWFMzYuBatUfGAiBjCQ6nH/5afsof9QIYRL4CAcmE3kUeAEA5EQAVw9wXGFrdUCi7P4arbyczpqnbmGTKdCxa7GvPcl0cjd7FXrFNrtR3UZLpYzlY2UTv0/XzvTtJeeoVZVcRjt2zaCiDHQfJcAJuSgJxKA8gcBIPPgryTyACKAHJOeDCH3APyWj3SLZYD4KquBaLKy9QybuLl2nazY4vDCRXjouo2GbRmMhPZi1uY1cJlOAf1LHFkgFRh2Uk+VUmUGqrL4Bd6onu3ivKVEzRjClZjxhYkLXjvIlbhbnTp+Nj3TNvZKlPMkX0I9RZybNg/aJuGLpBcpHLRq4IUfEyRRAxQXsoyc597rb9+VymtQdC8GIQFwpNp9V7EzWXNl1a1d8UZToUN6b2ISsuO6dQL9U3drlE4izRKd2cR9bXcxMtNGAGW1il0+iVeBolOr8PVajUYaOgq9B12NZwcJX4iMbpNIqKh46PSbM4xgyapkbM2jFFJs1QKCCZEycE6vYCPmYVhBooAHD4YqIHIVMFSARx7FEFDGDgPH3EnqCPAePVE7USAFTulDR7RJoD164crIi2jytyeRzPPWV9HgiQqGVEvqNygQVXBwOmmPS6LHt9k3PtrmcYdv8AiqTkB5TX0q1yxnbJ7a+V3EFGIg9CLZ1qkSQV00hky/yyTk9pq0lVoeyYifxFbl05m6NHDyFaSgQxu1S8I6q2htnDCWW0NQtg8gysyD2sYrxMXILTZuWlFF5xZvkbFddg5WNcz8xY2jCP/wBu9lhTStRTlX8YSzMUbK6Rd6nYZytsfvLOBQe7vQ5/tqtYIGljxPhWhZvfVCM2HWK5bsZazWbNlHsjJMjKAUkgq7jBEzYisshtrItOr1mVJWTOWbf9bdLMda6zFoyMSbu+Ys35CSSfZAy/lKyzlunnUw+UI9sSOOoiyycnD4Vo1lmBCaf4xxYnW6K2cM4tJvBgSGjCtp1yxg/FOdayNMzFj2oZEgzrnkI5tbq5E2ZtX7ADVyzSn66M2yeKwVlYIPXgQtkjAaS0OKp1I582UMHIZpV6vV6TWo2qU6HiKnVa7ERUNX4CvRTSKgYSHjmibSPjYevR6CDCJj2jVNNu2YM2qCCKJCETTKRMAC7oItE1XBiODKesi3M5bCYxkiAimHonI29ys/ULycxEyJisb6hKYQAQVe5HbTQRm2TjY/Ke+muxzJV+NrMKlQ2eeNf69EEEkE9krhc52rvsz0mJqzV3+Mrbbrda8wSU+jDGjIqa+ZTLhLfnEuYcX51rTW24evlbtsIs1iH8qpFuiDYGgSLIXUdF22Fdkb2Spynw3mdSIssfEzzJRBVo5YoHTcpJhoPZ5A2uHc3qNuExqph7eTG56HdJ+VAJslv2pxsSGSwVUIRI4vpyqqq4Vi8tSajdiyjqnJEiTuLA+TnUYhuu1wj0hUjmIU6xUlBIcSh48AAiAiAHAgGKXjj6PIB9vDkB6XB3O8Z2u0a6Fy3iuuuJzOurl9q+wOE3pHTIicHZoZdaq3Gadxki6QhLA2i8R2/IoGhrCk6Yj5FeMmqsyzjPHdzGmQqZmbHdHyrjeeJacc5KqsDf6TY27d2xYy1Xssa3l4F8mg/bspVJCSYPWzxNB4yQcIJiCaqCKgCmAYzsZnGq4BwnkPMNskm0TC0eCNIHWkCSItV5ORcN4aqRy/yto8kk0p60ScJDAu1aqmaqSKarr0G5FlU9d+31he34tw0pecvtHrfYLYeff5szW1nTR8pZ6XYr25c2Vphpa1tnb91Zalgwky7x3RFVZF0izrzJFFmgzQP6IQjnNo7213UxhrGgBnmHtYVGOwGyrqNcOnkbL5GmYwzfDWv9/rc/8PHWSi3ir2qeyymdFnNM4ixY0rvr/AvgaFBqxEVETpmDyEC/u0gMg3MqYREoiILCAn+sAExxMcBEQ5Hk3QXro6OjoDo6OjoDo6OjoDr5rJprJKoqkIqkqmdNVNQoHIomcolOQ5DAJTEOURKYpgEDAIgICA9fToH7D/8AX7f+/QKL18nf/CFtHmzT+TEzPGuTWFk2d1OaqA3SM8ikRMbOmHKFBxInhKXQcHEVxx+Fo2Ub11F+a+v/AJOnIA0kBaQhp3s5hnE2utl3nzdY3SCW+eZLZnHFkM2q1uyRnZliWyGaDjDEcpG0SEtVrlIPGS3z0i7eNCSxvRlLSmLV9HBOLi4sH9IW112Pz5pNENtN4V812Vgc3UheGvlYukBje41uhuSyxb7Gs77Mz9YkmtVm1Ea+ewVWLmPTsBY5kVzFP/gkQSj7TrQjub6not7dFxvbou+WpetsI2dvd0vG2SzaNMdETy6OMscMK+rjPBUPPKGQNP1LC8DUapKKR0aZzHuSxzEUQ3fQwVsnu8ghJ7drutf9fllBtVLwBhbJd6pOZ5J0/EposNgcn0aZjXMPLVNuiYidVxRfH1Ds3z2SJa2z/wCVw4IMYp1WrdBrcJRseQcPUK1U4xhGxdThItrFQcVCEKcjZCFaMW6UYiZMEz+LdqBUQMYx3AFMcgmX4RbvLESKiEP2wzplKJOFbLtcqJij/wDqGUqZjHH/ABMIj+o9fsHXeXDx4hu2AAEAQIAWLawAKA8cgUAqX0/YPtx0DSy/wl/yD/t170rf4/vOf8J7Yf8A1Jtd/KfR8f3nP+FdsP8A6k2u/lPoGcySyiDUyqRCnMUxfY4GMmAe/JlSkAVDJh/WKkU6g8h4lH360MzTpNCWq0OMr4AvDjUnP8jJSz6zZVxRTKN8RkqImxaqy8ZlWrSsWtVMiyCizNr8nulzipu31BNzL/hSVjRnZgHkcC+7zY/eJ7YQ/wCdk2tH/vU+qcFO8oUDf/Be2EcTCYTCpZNrVDCBuPInkpUzD4Dx/Bz4/wCHQZFQtzZqBsrzDe8OOEtbrMlFykVB5blLJAJa97ALxJUmdnmsZWsJx3I1FmkWQiVImv5fTo9xsPzRX8NwsoERLgw1s1RzfR9J8N7gYFylbG1JqeitsnbBjysvGb6ZTq2mltOohrY/fzddZyzubStEjB3Fsksq+kbwmVodSyItUzR51MxyHjXunZerqlKyxi3tLZGq71ZtIo1e9E2PtUKk+jCLkQlm8fMUZ+xaSzMr0wsZBBJN2zMooZoskYxh65+6toT3eMfd13EcLktHHmTdKbzbddlsxpVS81qw1BbEWHT5INRcVWlLMcsyznkxrUlLrNHnVLFE2FC3/Hxh5F3NqQzb4IOpHt44pv1PxTMZaznXi17Y/ZazP81ZmgnTtnYZDHctaR+OYYSibkm5k3kvjrFRnUjH0GMXmXrSFaysmRim2K7XBRgnVrZFAFDFIVNMqYmAxES+KXqiAAoBSiUogVMSlAnsACAjwI8dXToDo6OjoDo6OjoDo6OjoDo6OjoDo6OjoDqndt27xq5aO0EnLR03WbumzhMqyDhuumZJZBZI4GIqkqmYyaiZymKchjFMAgIh1UdeGHgph9vYBH3+32/P/DoFBa9XGM01z7sHpxNuFm1Omoe5bZamQfpopryuL0Fxe5txfQomLAa3SKDr5Lz+NahR4CTCufER1p5hGT1pHv1GWleFu41rvqdpTdN6coTit/ydu9k6R2VwzgSkLtLztrYcdZvVkLpgXA7+uoqv5dZbH1QVnY6KKq5NjWuokeM4CYRQfIkcZx/SGNadnM/aiY5T0tiCMtlY/YrH5WN8r1tg8eW+HxVKRVvNeIQLzNzVbfKUaWlEKo7tlDQl1WdlPFRRlIKT+WJ/D6d6Ydl3efVnMrraidhe29njaZZeztIvMuQZnZqGj6nVpuQYPYmr40xFVao0xPjGDp6LAkbSk6RVK+5rkO4dxcaVgzdOG6gSRXNdN7+9XLQuTtybDMazdrjIDgbVWNEI9RepZ6yFWIsgpUhfOVyrBEpWLgL1DTLuWtlJbXsikdORkSQkCgDY3pdAGuOuWFNTMQ0jXbAFEisZ40pEKgwga5DI+BgQYJNmx5KxzYCpKTtklABJWUsM2+fTky5TO6fPnS3qKdakA57ywFEpYXtflKIAAgWw7VgHBQEChwFS4+kBHx9vbn26+RTd5QphOEP2wxMb+MTWba83qe/IApzVB9QC/YgH5AgCIF45HoGpdHSt/j+85/wrth/9SbXfyn0fH95z/hXbD/6k2u/lPoGevBECJCUQAfXJ9xVABDxOIh+6ARH2AR4P+79vqEB8ekzbe9rWLvl5s20ujeRZfU3d6ZtVSvg3uDstrY4Xy3L05KYafK86YYiV3lFuKcylYnq8vcZKizFzB4AKlkfN67MrLZn3ebMAgMT2w+B//km138p9U4qd5PyA4QnbAKoUDeJyWHashygf+MCmCpAIAcQKJwEfExilEwCJQ4DW/Hvc4UjMh5T087lFKjNUM8v78rizHeQq3GX+Y1mzxX8hx9mTqU9Qcny7CRiq9IKxkOu7XgsnyVSkGqrxo3bRpRK+TR1T7PeaqtofrT3DNR8nT2WndW7XWc7fWWGT84M5B2nasW3WYsDHCJYFhHhIyjeuN0KmWPZoV6NRq7VCVTNBD8CQ6qe32xut3cr20xJP4J2Cxp2ur/j21N2yq0M+tO2TSSjnEYYRbz1fsMdWG8vWrSxFYSx1ig5FhMsU3DtJm+SSdOSqc4OFO0L3PcF95jFFhv8AJSOV9RLZlGn23LiLbZex5DxfJ48bRtiVgMaXCn53vx8p5UY0BR6k0ZBcarY25Tn+Ji13AmXWAOwnt34hvdCxHIZLzjAmgNlNjrJJ5rzrCOn7Wfd49s9uVPLpYRibcR1IvLJjzDy0lKV3HhXEtJtIiJcOko5RFJ4qCrBOrWyIBFDEKBABITlHwARICo8esRPzDzIkQxSgQpfEghxwHsHF06A6Ojo6A6Ojo6A6Ojo6A6Ojo6A6Ojo6A6pnrRq/Zu2L5ug7ZPWy7R40cpEXbOmrlI6Lhu4QVKZNZBdE50lUlCmIomYxDlEoiA1PXg88DxxzwPHP25/x/wAP16BPuu9yZaY5m2P0+sLt4pU0YO67a6sRSoIhITuIRcGl8wY6pMcwEa5TKLgGwWKiY/odekhroqQ8yiMOwdMWLtZpDem+zOJcEastNss02R7JTm/mSrdszi6tRFbuGRc6O8cZckHdywhiyzNajDWOfdRWKq1PI0laUdrL4uoTx8zi42dj46Sa+vgn9IZ1m2b2A1bxSXSqIWjNiIXYenEdZDg7fXKBYoXFj2KtSdygZK2TU7XpOVx9Kzv4WfWHH7J9IM7EMXHndwEh8uS9H4alaKdzPUhN1OVuI7cFwydOwrKMsmQbdd9uHgtwBNM8wwxjSwry1KwTSp1+k3kHmMsNRFOx8gqwh0WsARvCRRGgbjtddNi90Csp7dd2vhTDqCgWmj61YQyXdqhelX0mcrtiGxuSabMxjp3ZKiwO6rslQ8eXeew5ajTMq6mGssaJrzhJj9Xr8HToOFp9EhYWp1qqx8ZDR1Sg4lpDQUHBtmxUo1hW46ObN4po1Zt0km7diwTSYNW4CmVNMxEiAvMi3eWIkVH5P2wjplICfClk2uVExQ448xUqZhUH2DkxxER9xEREevoDvvMB48Q3bADwL4l4sW1geJfb2DipewewewcB0DS+jpW/x/ec/wCFdsP/AKk2u/lPo+P7zn/Cu2H/ANSbXfyn0DN5JZRBAqiSYKGBQoD5gIplKI8CZQpAFQxRHgoAmUxvMxBEAIBhBf8AlvSNs6ta+VdX77/4R82SD6Zf3aw43pVIJWsysJaQLOO47LdPeRDisWCafTraNFTKjuGe5Phos0zHQNgbtZ6XbPcHF93mzBwMT2wh/wA7JtaPH+XNT6pgU7yhSmD5N2wjCYTCJj2Ta1Q4+RvISiZSpmESAbgSkMIkDgOAAADoL/TNvkrE7s+BN38ZG1auUtXJqrs5m0WqvoYbz+xVElStllw7dGVhfmiYiQkJVEarT8kPKpliUh5Ikk2qqpYmbWjtbNT8/o6iap7N4pzHYWtck9CLzZqLV6ZIs1niFc19nrG9hdHYKTmKu0kDSv4xoJac2CSTfyFnSIudzeFmMmDw5MpyjibuiZsrhqXl7FfaayHBGUCVj4m3G2TniQj9Bk5jiWCvqSNIchA2ZkjILhGTkSozlI46xzsnaAiY3SEcRaMd4Kh92zGNfzotRssaUWHJ2JLtld5WbpUJamTtUxThy6U7ElKtMXkmVY5tyWXHbKcj4+VeXOBnVLXZIxjc5J1NzDFrNJB1L9vjCtxxThhW2ZfZOGuwWe7BK5uzg2lnbSem6jccgPHNnVw6ncUnEi7sFJwqtMvaFjlF1MSKEPV2CDKPFFsPgO+nVtZ8eofx8QAplPICAIEMsY3kuIAYAN9KnJQEQ+oB55Hq5dAdHR0dAdHR0dB//9k=)![ref3]![ref10]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEABcDASIAAhEBAxEB/8QAFwABAAMAAAAAAAAAAAAAAAAAAAYHCv/EACcQAAAEBAUFAQEAAAAAAAAAAAECBQYDBBESAAcICRQTFRYiQRgk/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/ANb2YWgQcyXm53wOtXcAYYu1Rnlw7Sy21HeLMhvHi2mFNayB4bPdoSCX2w5LlzFpSlDqjSowpK2yAip0lFHcG3OyDEloRrIeq+yGWpAG0hQYXqUPgVGgfcMMBY2XO30XLp5ozy/bO4M+xRu40a2Y+pLyhmqYqCVPJVVdDFmSXNGSCeFQT/6oXGVJSSm/fj9MzDDAf//Z)![ref11]![ref12]![ref13]![ref12]![ref14]![ref14]![ref14]![ref2]

<a name="br18"></a> 

18

Fast Approximations and Coresets for (k, ℓ)-Median under Dynamic Time Warping

and further observe that for

dtw | (τ , c<sup>opt</sup>) ≤ dtw | (τ , τ) + dtw | (τ, c<sup>opt</sup>)

∗

∗

p

X

p

X

p

X

i

i

≤ dtw<sub>p</sub>(τ , τ

) + dtw (τ, c<sup>opt</sup>)

∗

p

i

≤ (1 + ε) dtw<sub>p</sub>(c , τ

opt ) + dtw ( <sup>opt</sup>)

p τ, c<sub>i</sub>

i

≤ (2 + ε) dtw<sub>p</sub>(c<sub>i</sub> , τ .

opt

)

[P](#br15)[ ](#br15)P

And thus it holds that

conjunction with Lemma [22,](#br15)[ ](#br15)Observation [24](#br17)[ ](#br17)and Lemma [25](#br17)[ ](#br17)this yields

dtw | (τ , π ) ≤ (4 + 2ε)∆ (refer to Figure [6).](#br19)[ ](#br19)In

∗

∗

i

∗

p

X

i

<sup>τ∈V</sup>i

X

X

∆ = dtw<sub>p</sub>(τ, c(τ)) ≤ dtw (τ, c(τ ))

∗

p

τ∈X

τ∈X

(

τ, c τ

X

≤ (2m)

1/p dtw<sub>p |X</sub>(

∗))

τ∈X

ꢂ

ꢃ

X

≤ (2m)<sup>1/p</sup>

≤ (2m)<sup>1/p</sup>

≤ (2m)<sup>1/p</sup>

≤ (2m)<sup>1/p</sup>

dtw<sub>p |X</sub>

(

τ, τ

∗) + dtw<sub>p |X</sub>(

∗

τ , c τ

(

∗))

τ∈X









ꢂ

!

X

X

dtw ( ∗) + dtw<sub>p |X∗</sub> (

∗

∗))

τ , c τ

(

p τ, τ

τ∈X

τ∈X

!

X X

X X

dtw ( ∗) +

α

α

dtw<sub>p |X∗</sub>

(

∗)

τ<sup>∗</sup>, π<sub>i</sub>

p τ, τ

i

<sup>τ∈V</sup>i

i

<sup>τ∈V</sup>i

!

X X

X X

dtw ( ∗) +

dtw (

p τ , π

∗)

i

∗

p τ, τ

i

<sup>τ∈V</sup>i

i

<sup>τ∈V</sup>i

!

X X

X X

≤ (2m)<sup>1/p</sup>

(1 + ) dtw ( <sup>opt</sup>) + (2 )<sup>1</sup><sub>/p</sub>

p τ, c<sub>i</sub>

ε

α

ℓ

dtw<sub>p |X</sub>

(

∗)

τ<sup>∗</sup>, π<sub>i</sub>

i

<sup>τ∈V</sup>i

i

<sup>τ∈V</sup>i

ꢃ

≤ (2m)<sup>1/p</sup>

(1 + )∆∗ + (2 )<sup>1</sup><sub>/p</sub>(4 + 2 )∆∗ (4 )<sup>1</sup><sub>/p</sub>((4 + 2 ) + 1 + )∆∗

ε

α

ℓ

ε

≤

mℓ

ε α

ε

.

◀

▶ Lemma 27. Let X ⊂ X<sup>d</sup> be a set of n curves. The metric closure dtw | for all pairs of

p

X

ℓ

curves in X can be computed in O(n<sup>2</sup>ℓ<sup>2</sup>d + n<sup>3</sup>) time.

Proof. First compute the value dtw (σ, τ) for all pairs of curves σ, τ ∈ X. This takes O(n<sup>2</sup>ℓ<sup>2</sup>)

p

ꢀ ꢁ

time. From this deﬁne the complete graph G(X) = (X, ) on X, where the edge weights

X

2

correspond to the computed values. Clearly the metric closure of dtw corresponds to the

p

ꢀ ꢁ

weights of a shortest path in G(X). All these values can be computed in O(n<sup>3</sup>) time by

n

2

n applications of Dijkstra’s algorithm.

◀

▶ Theorem 28 ([[21](#br25)]). Given a set P of n points in a metric space, for 0 < ε < 1, one can

compute a (10 + ε)-approximate k-median clustering of P in O(nk + k<sup>7</sup>ε <sup>5</sup> log<sup>5</sup> n) time, with

−

constant probability of success.

▶ Theorem 29. Let X be a set of curves of complexity at most m. Let k and ℓ be given.

Let X = {τ | τ ∈ X} be a set of (1 + ε)-approximate optimal ℓ-simpliﬁcations. There is

∗

∗

an algorithm with input X , which computes a (10 + ε, 1)-approximation to the k-median

∗

problem of X<sub>∗</sub> in (X , dtw |<sub>X∗</sub> ) in O(n<sup>2</sup>ℓ<sup>2</sup>d + n<sup>3</sup> + nk + k<sup>7</sup>ε<sub>−</sub><sup>5</sup> log<sup>5</sup> n) time.

∗

p

Proof. This is a direct consequence of Lemma [27](#br18)[ ](#br18)and Theorem [28.](#br18)

◀

We next show how to combine our ideas with Indyk’s sampling technique for bicriteria

k-median approximation [\[29\]](#br25)[ ](#br25)to achieve linear dependence on n.

![ref14]![ref14]![ref14]![ref11]![ref13]![ref14]![ref13]![ref14]![ref14]![ref13]![ref12]![ref13]![ref12]![ref14]![ref11]![ref12]

<a name="br19"></a> 

Conradi, Kolbe, Psarros and Rohde

19

V<sub>i</sub>

V<sub>j</sub>

π

∗

i

τ

∗

c

opt

j

π<sub>i</sub>

opt

i

c

Figure 6 Illustration to Proof of Lemma [26:](#br17)[ ](#br17)Assigning τ (the (1 + ε)-simpliﬁcation of τ which

∗

lies inside the Voronoi cell V of c<sup>opt</sup>) to π (the (1 + ε)-simpliﬁcation of the closest element π in V

i

∗

i

i

i

i

to c<sup>opt</sup>) under dtw | is at most 4 + 2ε times as bad as assigning τ to c<sup>opt</sup> under dtw .

p

X∗

p

i

i

5\.2 Linear time algorithm

With Theorem [29](#br18)[ ](#br18)(and by extension also Corollary [45)](#br27)[ ](#br27)we have ran into the following

predicament: We would like to apply linear time algorithms to the metric closure of dtw .

p

However, constructing the metric closure takes cubic time. We circumvent this by applying

the following algorithm, which reduces a k-median instance with n points to two k-median

√

instances with O( n) points, simply by sampling. More precisely, we will apply this technique

twice, so that we will compute the metric closure only on sampled subsets of size O(n<sup>1</sup><sub>/</sub><sup>4</sup>).

In this section we want to analyse the problem of computing a k-median of a set X in the

metric space (X, ϕ), where ϕ is a distance function on X with the guarantee that there is a

constant ζ such that for any x, y ∈ X it holds that ϕ(x, y) ≤ ζϕ(x, y), with a linear running

time, and more precisely, only a linear number of calls to the distance function ϕ, and no

calls to ϕ. By Lemma [25](#br17)[ ](#br17)the results in this section translate directly to ϕ = dtw | with

p

X

ζ = (m + ℓ)

<sup>1</sup><sub>/p</sub>.

Observe, that similar to Theorem [29,](#br18)[ ](#br18)the following lemma holds.

▶ Lemma 30. Let X be a set of n points, equipped with a distance function ϕ that can be

computed in time T . There is a (10+ε, 1)-approximate algorithm for k-median of X in (X, ϕ)

ϕ

that has constant probability of success and has running time O(n<sup>2</sup>T +n<sup>3</sup> +nk+k<sup>7</sup>ε <sup>5</sup> log<sup>5</sup> n).

−

ϕ

▶ Lemma 31. Let X be a set of n points, equipped with a distance function ϕ, such that

ϕ ≤ ζϕ for some ζ > 0, and Y ⊂ X. A (α, β)-approximation for the k-median problem for

Y in (Y, ϕ| ) is a (αζ, β)-approximation for the k-median problem for Y in (Y, ϕ| ).

Y

Y

Proof. Let C ⊂ Y be a (α, β)-approximation for the k-median problem for Y in (Y, ϕ| ),

Y

and let C<sup>opt</sup> = {c<sup>opt</sup>, . . . , c<sup>opt</sup>} ⊂ Y be an optimal solution for the k-median problem for

1

k

Y in (Y, ϕ| ) with cost ∆

<sup>opt</sup>. For any

τ ∈ Y

let ( ) be the closest element among

c

τ

Y

Y

C under ϕ| , and let c

<sup>opt</sup>( ) be the closest element among

X

τ

C<sup>opt under</sup>

ϕ

\=

ϕ|<sub>Y</sub>

. Then

Y

P

∆<sup>opt</sup> =

ϕ(τ, c<sup>opt</sup>(τ)). Overall by Observation [24](#br17)[ ](#br17)we see that

τ∈Y

X

X

X

X

ϕ(τ, c<sub>Y</sub> (τ)) ≤ ϕ| (τ, c (τ)) ≤ α ϕ| (τ, c (τ))

Y

Y

Y

X

τ∈Y

τ∈Y

τ∈Y

X

X

≤ α ϕ(τ, c<sub>X</sub>(τ)) ≤ α ζϕ(τ, c<sub>X</sub>(τ)) = ζα∆

opt

.

◀

τ∈Y

τ∈Y

▶ Theorem 32 ([[29](#br25)]). Let A be a (α, β)-approximate algorithm for k-median in metric spaces

with constant success probability. Then for any ε > 0 the k-Routine in Algorithm [1](#br20)[ ](#br20)provided

![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCACdAaMDASIAAhEBAxEB/8QAHgABAAAGAwEAAAAAAAAAAAAAAAEFBgcICQIDBAr/xABSEAABBAEDAwEFBAUKAgUICwACAQMEBQYABxEIEiETCRQiMVEVFiNBGDJYYZEXGVNxgZKY0dTWJJUlM0JVlig4Vll3scHYQ1JiY4KHobjV8PH/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A+/jTTTQNNNNA0000DTTTQNNNeKQ680S9naff2o2Kp2oCpz3kZ/FzzynCdqIn7+dB7dNSspxB6huIaNCJcqjfINq0qIaq53fGjqmKtL2jygH4/LXBmRLQV9U23+DcJXGG0RUEeEGMjfcvc4nPKO9yd3C8Amgm+mpKNhK9B1w44tu9huRo7ziNG40HHBur2L6JqK8uN8GjSoiKR93KS1nMsbcDuLI8f4VV7CC5riE2+fhNfx07VJPKh57flyugqzTVMffDGP8A0kov+b1v+o0XMMZ4XjI6FV48ItvWonP5cr7x4T9+gqfTUggXTVjHdkRJECU13oyy/AmszWHJIoXvEZHGeW1OOYiBL3ckhKSiCpwvczMmOcuOtiwTDZBJhiPrGrxKPYbT/IdzAoLnYXpor4r38N9vChOdNS5JTqtLwokaIqepx2ivH6r3Z54aNEVeVJeE4Tzz46RnOOqIsk24iMA6840nqEBOKPpALXwoSPD6hiXcnaja/CqL4Cb6aaaBpppoGmmmgaaaaBppqWyZTsdxSXtRgTbF03E9MWgMV/FE+S9X8TtbUO0OFNPiXjhQmWmpYkp0vUETaXudNtpxPICoKSEBfLlxOP1fCeCVF4ThfJY3cSnaGVbWECtioPpk9PkMQ2TkKqKnD77gAnCIai2iERonPhBXQT7TWsrNva69COJ5JlW3VFvpje7e9WL3NnjQbB7NmmcbzZHkdLLdiXNFjOHida3a29MEWfPsoYWbYsV1ZYyUfc92QHLCZF7T3qZztWpvSz0FZhfUFSHuuWT+rLOHela0W0eTvhHhVO1hm633spjaFz3u0WbU+6SlYi+7Pet6wBugsXPTYRTESZIkF7lOSES8AYD/ANohdVteOR4TkufGvGrpNtmZuPAEdG0kuK0gE+jbfxPeqrvHpcCpmS+BDnlePOtBNrlftX83WTLm9W3TjtrTZUpTXts6/ptmXGX4dU3CK7b4dj+5/wDK1WDeW2JV0iVVV24q4jWfadhAj3f3YgpJ9wYpxzo+yeQwtfbddftBsloDaZSRjF71JFYU2R0c2OoPU1tXJh7PrVc+uechWkFHvx4jr7COj3qSB9Bn3txL/wBJMc/5zXf6rWF+8/tPfZ99POe2O1+9XVxsxttn9TErJ9lieRZW3GtIcO5hNWNZIeaYYkNi3NgvsyWVR1VJtwVVBVeNal/5qv2dn7Iezn/IX/8AW6yk2u2D2Y2Uw+Ft/tVtriOD4ZWyZ8yDjtLURm6+NKtJTk6weaGQL7qHKluuPu8uKiuGqoiJ40F8P56z2Un7enT5/wCMy/0Wn89Z7KT9vTp8/wDGZf6LVJfdvH/+46f/AJXA/wBPp928f/7jp/8AlcD/AE+gq3+es9lJ+3p0+f8AjMv9Fr21ftkfZc5Jb1WP0fXNsFZXN5ZQaengRMu9SXNtrSS3Ar4jIFEEXXJUqQ0wAEQIqmgqSJ51Qv3bx/8A7jp/+VwP9Prx2eGYncVlhTWWOUkurtoEyssoTtXC9GbX2EdyJNivoLAqrUiM86y4gkKqBkiKi+dBtMXMMdUG2iySkUwUmicWbXIwSiBEJi37+vaPeI9pd69y8BxyXKTeNYMTmAkQJ7E9HwQUkQzadbQmxRHSJxl10EQXhVFFDLtXiOvlVNPm/L2VPs6iRUHpD2cb57UVW6B4eUEkPhf+MVOFUeC8eRVU1VNf0O49hsQca2K6gurDpm22YUjptmenfeU9v9qqCZLc95tptDii4/arCnZDZuSbe7d9/c9/tZ0yT2tet2CH0NNEpC0rSP8AarxCBstK223IFDSS5KTvJBB10XO1fPBGI+eedTiN6XpfhCgp3n3IicIrvcvql/8AiPuXn8+edfPxUZd7Unaswn0/UFst1KYZhiuV9FsbmW0cva7cDcerryKvp6PJ+o1M7zYIGYVFcLF/kOWJtdPTKrCknxkp6lLxH4F6sW9qduzt3ENer/os3N2+GdLjyoOVdOls51K7aYxjCN8Xd/unmruP7YFgi0LynLnV4092LNHGmWPvarG93cDdLprFDp/61+l7qmor7Iunnffbjd2oxW3Yor6xxG8jyYka5mQhsItYRvKwfrFGcbdUgbP4UL4VXhNZLsWJSmW346tvCgEbiNL6iOEK9jjTJIg8o05yPqcJ3ECj2p3coE401K0lvoiJ2dyiSK6pj6Rdjid4AAfH3E2hCDi8p8QlzxqYNd/YiuGJkqkqKI9nAqqqIqPJcEIqiEvPlU54TnhA7NNNNA0000DTTTQNNNNA0005RPmvGg6nnm2AVx0kAE4RVX6r8k1wWSCKSKJpwgKiqiIhqfPAgvPxEnC9yeOPHledcZakrKqBfqqhGggjhmCIvIgncPBKqovdyvCIvjzrXR1He0n6aumfLJm1lrZ3u6m/VZAgWIbD7L0beYbuBj1wMj0swn4sM+ElXhEWVFjR7a+Ka+cMpkRVgn39ug2MHKaBon1UlbEVLkUReRTjyKc+UXlONYgdSnXP0kdKkzH6vqI33252tn5kVkOOU+W3jMewyZKFYrdtBqITaOG5OYkWNbHRhwmiN6WyHIqvOtM+YbwdfvVLaMS9x9xj6N9t5TSHC2y6csnk2+61ZkFVyMa4teok4dAErFLpmY+Vvt//ACZMKxKiQD+8T3o8FINt+n7bDbKyyHIaask3mdZ8DFruXu3mQBle4We5JBJxW8kzLIZAxftTK7ZyfOflyW4sX3twnjURUURAyZzT2tG8uZSo73Sn0Mbi5HGrWkTMy6qckPpceZWX5oFwRGcX3VTK2H2mLIraWS1P2WoVg+jLSw5j4+326/tMd4DnDl/UVtxsft1nBK5cbabX7aSrPdzbOql8uLQ4f1KrmlWNnkVIi+6xM6Xa2uSWhOv/AHfjo56Q3cAG2u5ttphkVdcUWYrnbGZ7eEIQjdv/AAzpKqetE9R33VRAUedQ+R7NBjIWw27RgrD3X17QObGMSCRFn9RbkiHNiGnYdfJjfdMFJp1lXGZLnqf8Q04o9jfz1bU/ZodBneSr0vbUSVUlIn3aF8XHSVeScNPfj+Iy5Jfi+arrOXTQYMfzaHQZ+yvtJ/yN/wD1en82h0Gfsr7Sf8jf/wBXrOfTQYoUXShF23gN4zsD1AdTnTRtuy45Mj7UdPW67m3O3sW3mr3213Fo1prlWLa1eBt62l+8H7+622fpsdnC1vRsdeO0Quwtieu3LbiDfKknK3+rvD5HU7OamwUIKv7ivN5htn9zoD4SpxXcMhuPtN8Kx5Ho3uHY/fjTQUpRdePtM9ta1rEco2B2A6q56EUr+Vmq3XldMaPBMT1Bx1ra8sJ3eXuqlbRs7hcxX7XTl73CD2emWRGI+2V2xqWKGs6mNiOoHp7ua5hhrdPci4wFLzpl2+tm2eyzumt6isqc7nAWZ3/AVOWFh1ctkk2E79lRvePTCz7iCqJ3GAcF3CJve7o8SIvDKP8AY56RH80P0z/V47V55TqcYCfHUJUJoWHnY7suHJEZrLqzG3nS7JbgAM/hW15rSjMoP/Xe8fgdhhsw2N9pT0GdS2XycC2H6qtntz8xh0cvJZWP4xlDEmwZoIMmFDl2ptvtxw90YlWMFkzQ1LvktcCqFymaMO0r7GME2tmRrGI4qoEiA+zKZJR8EiOsmYKor4JO7lF8Kia+Z/dDYDZrevHWcW3S2uxfNKCLdxraNR29cyLLMqtYm18K19WOLThmxEnyWe/hBFJKiokpIQ22oukxvbeAeLdP3UZ1OdMe3Ud552r2v2F3cd2923hynzU7S4iUK0luqW9zIQJVlM964ff7i9Ie/wAB9V7kpltTFVVSBQTtFO4iVxFUUFOfK8Iqqnjwi64JNaU+xQeEu4kFCBEUwFeFdD4via54Tv8A/tD486+ZKryH2jm2deOLbL9b9HbYEwrsqK91QbPSt+d0pdk+YrchJ3Gb3GwUZ9d76oFjlR93WfsuCJxVlyuUeS4OO9aftQdnWZDeVYVsB1rSrJxh2JctTJfR9KwBYYG1KopFCsTfb72uXLj4SnbZLKiSlOr909xsPtL14gfRvyn1T+Kacp9U/imvnx/nQvaP/wDq7Niv8bEr/wCXnT+dC9o//wCrs2K/xsSv/l50H0GEYjxyvzXhETyqr9E4/dyv9muBvgBo3wREoKaIKIvgSEV/NPPJIvH05X8tfOnk3Xd7ULdyC1juMbRdPfRxYRJI2p7kWOey+qtm7jstuxzxAdvUxjZda1yaUkbMck+9Ev3T7JWD9kv/AGj7zFtraTfaHbsQhx/fLrfKpxFh4JlS50s7cyemncI8lbE22RucyczrchLXFSr3bJJlAlPEWXarU2P2g19nrHkB9ME3JsfrXW48+5rIch0HXGWJc2NGdeBhwWnjaCQ62rgtOGIGo8oKkiL89awN6fbF9Du1Nzl2OUOd3u/24e3ttb1Oe7XdOuON7l7i4REx1JjWUXmS0aWtO3CxrHrKG3WXNsMx0YVlLr4yMupI9QNS170WbSbpPMS+pXItzOr6bTKUDFHup/LXdy7fDIElsztYOFvjEpfRYfkMsvWzCgnvkuMxM72vQ9I7/QdsNv6TA3Nucaxyqx3DGcJfwikg0MJmskJRnAbpYrRmiPONOM0hPNtd5v8AqSPSfLyHYoVpD63uv3q1xinyzpyxjZrpc2ozqLD3E2x34zazk7/5Jnm2dtHGTh1bd7ENQdrwwW9v6Ozr722sB3IyMcfmQHqAYtmlgllFpm56L6jeAFe6vt1t1erMbOSl5Z7fbm5O47sPjucPturNzHbbalWpJYcsZJNnVYvUlk9t93cbuJlOs6x596WhfZgZRaW3TC9hUgUYodid4d5OnPBmWk9Entvdi84s9u8Gcuk5NLXKWcdpYbeRXY+6BOszkPjBjo/6YbEtBR2I4BhuC0OOYziuM0lFQYjWV2OYzV1lbFjBS43VQUhQK2G+gk8MWO1HiNCwpdxNtohPKqL3Vf8ArvNm82Lzgt8tTHfjOC83wAkwz4QwJj1GWw7x9AHETlzs+LlpoAuuihECPMrIMieZYkI2AE0ittOAStkpLIaVVkpwPLhkXP5a85NL3AhIj8dt5uS0yn4DrbhMGDjRSPxO9pk3FRtEaFDEETgedejTQNNNNA0000DTTTQNQXyip58ovyXhf7F88L9F48L51HTQcCQ1BEVQeFfiKPIBDiEQ/AHrNISK+6IfiNv9zatyRbe7C7exeowUu4lbcdJ1G47yuyEJ44fAm416xNKKgpALRgrS+qKk4vCfAvo00GMu8HR503b6XFZlm420uLW+4OMxPdME3BixArchwEhnDZwrLGjjj2xLuqniE+JZl6i/aLIy0ZFC9FJdUxet/pmkjcbB7/ZX1N4iyrdjlOz3VvmL95lmdXDzSVrcfDN+24ips1jNMy+N3LpV24zhb+XVuRPfa37W94g5V64mnIknah8iqdpfqlyi/CXz8L8l8fJdBdzpq65cA32l2OCZPhmZbFb2YfBiWVxt/u7VtUC3TJz4dHeZRtZavSkLPduhyOa5SUmWpAp3LRXoTv2RF94VtvPaPJaFryvHJOKvxd/by6oIJL44LuVB7ePhVe3leOV0gb0dO2zm/sKmibo4izksjGrNi/xO0WP3ZBgeTQ61xqtybELhDE6+9pXSSxpn1juBFkMsuq24oIqzvZPqn3o6W8jq8G6yM9rdztms0nV1Dt91E12BO4bYbdWJymKTD8T33B7JMiDInstaGFFk7ypPpUttybiuxz7lRxt1s4obr2ngdTkELj4k5VOPIl2qnzXyiouu3UthyPVcQm+xyK/HCQ0+BCqcOdith3CpAaG2SOiYkqdvCJynxLMeU+qfxTQR0000DTTTQNNdTykIoQkgoJIR8j3KQJzyApynBF44LzxwvhdS9Zr6Ciq2okqk4gtj6yq0HHLZJ3B2Ol3J2/rc8F48aCa6tbu5ufgWy2DZPulufllZhuCYbUv3WQ5DePizVU8OIBEU18zIUAQVe1UHvNwzBBHwuseOqrrk2q6SomPR8wiZrnec5olj9zNp9qsYby/cy/hUiREvcog40lnANzEqF2zqG763KQCVv2nCIWJPrKgaaMgqN2upPNoW8nVDl+SPGlzDybEumSlyxy02B2xfqCcdw6adO7XRhzLczFwmWLUjcYmaELNJnp/dqH7v3OhcXdfrM6m+racFRsVPzPo82BZaKl3Evc+wt6q6lcwk25IdhM28Ju7jNbSRMfCBHXHMzIsy++A3sv8A6HpPsBPtKl9t9rcf23qnYcAshyHIbm1scozbLcwsQyHO86yG2Jp29yTObcosR2bbX8gGpNhO7QBH2lVI5d/w3LXyhB8mVQBRhOEZBpru9JgAT5MtdxemCqvbyvnUDRHT9R5PWMmn2nSd+JXxkq2Tyvr47yJWgXnx5T5aAn5eOP8A4fXThPlwnC/NOPnx8tR+WmgcJ9E0000DTTTQNNNNA0000DhF+ac6f/3+Hy/hppoGioi/NEX+vTTQQ4T6J4+Xj5aivnnnzyvK8/mv1X6rppoGmmmgaaaaBwn0Tx8v3c/P+OuowMQNwVQzdeYjNMqvAIyfJvkrnBei+itp6K9h+FUvy4Xt11maN/jm4rARB9U3wT1T9IzBk2iLkPdmjJ0VVrh33o0BFNr56DHn2e/pbfZn1hdN4E5YLtpvrP3abyJzglso/VHNyTdpmnVrlfc5GJg6lNMl+tITJH0W193qez3JdmWtXfTrYRdt+v3qbwTIm1ZuOovbza7d3bFqC13MfcrZKgr9t80dvmVIVx+zeyPNaM6mqErBJUBZb5S2ljI27tE0DTTTQNNNNA0000DTTTQNNNNA0000DTTTQNNNNA1I7+lqslp7PGMhiRpmO5FCmVV9FkNjIhT6KZDdjT6udFJEF73sSJsvjbVmMZPD3E2glPNQVELtRURRFVPsVEUCdVFEXiFfm82K9oH47UREVF40FmunzfnIukzdrbHpF3EyG1yvZjfK6m4f0p5LNdO4zDbO5oMWtc6mbUZp6px2pm3kHEMayB7brM0lQCxuvi4pt8mPWhB9vnuoYGQTbSiKg0hJ2F6vqq42K+HngQQ9U30/E7u4e1TQvi486aN69rqbeTbPKsBvEbArmuIqa9drkt52M5fBJuZhuY10JJEJyRcYrkcapuatG5kQhkV7CK8IovOTnszeoi/6luknAstzBuamZ4NkuebFZXJuLZZt1ll9sJmd/tFZ7g3ilFjuxrDciXhj2evVDgPfZL12sAbO1GMlhJDYzpppoGuhyS2z3+p3AIIK96p8Bd3dwILz8RJ2/EnCccp58679Smc6iH6aGAOkbLUc+z1vSef7+0nWu4FVv4F896c+fpoOyZNYbbVw1UQbUHPUVUEFRUL9RfPeQ+O4OB/WHz51g11gdYFD0041U4zQUsTP+orcMbFdr9oolglc9fyq0YoW2Y5POCNYPY7txjDtlVfejJ1rbD3H7TrxCA+r5EHX1pdXGP8ATRRUmN0tB/KXvxuW1bjs5tFBsEqpOWyqlIQZDkeS2nuti5j+3GMu2dIuXZG3WWv2Z9p1gFBd96Qg1QYPgl5EyG+3R3VyVNxd+c+KA5uDuSkEq9ZUCuWU5TYlS1DsyxPD8Cxw7CyHGsbSytUrVnWBHOle8ILQduE4bkLN7e7q7pZIueb7Z4Fc5uHn7cVaxuzagFMdr8Nxeudk2J4ttljpz5yY5hwzrJYHv80ztJHroLdzE5XkzQPVdJXXyEe1XHi/XccXle9wuE7j8d3CfCnGo6aBpppoGmmmgaf/AK6agvlF8onKfNV4Tz+ar+Sfv0HWTzQCRm4ANNm40884YNssOs8eo0+6ZCDZhz8SKvCcL510w58KyjMza2VHsYcgBcZlQXmpMcxJOeEdaMh7kRU7h55HlOdauPaC+z+3g6zkZiYJ1j7pbA4dY4nMw3cDbmDXPZht3mNaZwyg21fTtX2LrSttIxKG/tlkWC2hy64vd4vu/DtoOiP2VOa9P2zdhttnXVf1MV06v3GzmdSObPbsPYXiVvis+xaKgvHsXeqLw6a9s4jav3EFLacLD/a2MhxB7lDdfyX9G5/dT/PTkv6Nz+6n+esD/wBBUf2weuz/ABAr/tPT9BUf2weuz/ECv+09BnhyX9G5/dT/AD05L+jc/up/nrA/9BUf2weuz/ECv+09P0FR/bB67P8AECv+09BnhyX9G5/dT/PTkv6Nz+6n+esD/wBBUf2weuz/ABAr/tPT9BUf2weuz/ECv+09BnhyX9G5/dT/AD05L+jc/up/nrA/9BUf2weuz/ECv+09P0FR/bB67P8AECv+09BnhyX9G5/dT/PTkv6Nz+6n+esD/wBBUf2weuz/ABAr/tPT9BUf2weuz/ECv+09BnhyX9G5/dT/AD05L+jc/up/nrA/9BUf2weuz/ECv+09P0FR/bB67P8AECv+09BnhyX9G5/dT/PTkv6Nz+6n+esD/wBBUf2weuz/ABAr/tPT9BUf2weuz/ECv+09BnhyX9G5/dT/AD1BREnGT5dZ7AeR5WvgfLwpt8uJyhMoop3wlT8c/Te9Zv0ew8EP0FR/bB67P8QK/wC09cy6HSEAZ/TF68yjgByzL+XjujQ32jFtlHXPuuncj4umq/CPanj4vnoKhlTY2He0t6ec3yYmaTGM76Yd19pcZvZIe6wso3Wv9wNs8sqMIgvKRJNyqfjWKZRetVqembFdR2pIbqMKWtt70hmMk8pDgRwrO5Zzr5g01HbB0WCecdcIQFoXzbZUyJPxHAHj4vHzbdR/R0tJu90NV8nqo6u5YZV1bU9TPs7ner7TPEFl7PbuWbOT4tLXG4wY9kKvQQqW7c0kCEC2nwvRIpaEObvUf7KqVvpshuDtAz1vdarEvN6aPVpF3K3dc3BwdixamxJqu3GIt1ONHlcSQMR5xmAlzULGkehIWQ97v6bgbY4VrXWCuLDmxZIMOEzKdYkMvtQ3UVUbblGy44jJyEEjjiXlwBUvHGveiqvlRIfK8dycdw/9kx8ryBpwQL+Yqi8JzrQN7O/2Le4Xs/8AI6+/o+uvcrK6yRuZRZXmu31fjkjG9t9xcerMJy/HZeGZNTO5dee/5HAvr2puYlwr4jWsUc6ItdJWUkmPv4EEbTsExNE+XY16QDz57A+M/VAEXtF/4PXREd9Nvu7UDlpppoGmmmgaaaaBpppoGmmmgaiqEKCqiXxoqt+PLgivBkCc/EjfC96/kgr9NQ1wT0vVaZkK/wCgcn3knIpczGTbiusNCwKp+EhukCm/yXbyqemvPcgcCkMobLQOtuvSnPSiMtmCuy3EFTcGMBECvKyAuG8g+QBpwlTgV12koiTqKqKLHHrOCvLbfJo0qGX5KL5IwXjw6vb+/WDPVF0YTepbeLpk3ZPe/dzbtOnzN7DJ5+CYfeFjdTmIyaS/pooigpJGDZ+lbtSLnIjZmpd0TVnjH2bCW3S1hZyuNIBso0TDTbBkUlxtj3ZHy9ImWnvUVx30O4VAHnOHfWkkTSoHqd4hy0000HFztACkOOo23G7XD7j9IC5cBsEN/g/RETITU/Tc5Qe3tTnlLGdH8uw2z9pfvXtnTurDw7qF6cajqUy47Vv/AKQl7mYNmGIbKUsbFn+9pmvx1zB2Pe5leDMqRb2AleE5HF42AvqitibZOk622JIZOtp8LXYqGhuqqpw2iinI/wD0n/V8j3dyYP8AUdNDEeq32bu82Syjqtudt+qfK5O6G4Ermux7FKrPNhd0dssCfy6WhuDHjZFn+X4fh+Px3SJp+5u6dkeHXAXQfSBpppoGrQbz7rYLsdgOYbq7k5JU4nh2GUj99cXdtI90hxo9cy6aBJeQXXHFfcIGmWWWXXSVTUGy7VTVx0mvg44D5stKhctAq+XF/wC00nCeBFe3td/WPuXkB7fi+eXrI3jl9YnVBK2KjpUtdPnSBneC53k2R1D5zLDNupbH3r5pcFdmkxFGnDY5sCXM8f4sm8p+/wBTKcipWoH3wLY4TdbkdSGTY91adRtJNx/da9x2xZwHaq6cIo/TlhWTOMPWOB4u0QqM26uAr6ss+zMGaV7Lhq8bSTSRPslrvvxwn0+fz/fogqCGiRhhmpI3IaDhEdFlFWOaoiJygobvav59y8J9WgaaaaBpppoGmmmgtbu7vbtJsHjUHNt6dw8e2yw2ZkNXjIZJkc46+Kd1cDJWtq4b4NPKVlOWI+kVpUETFp5VNO3hblNvK9HadMkjpw67IhxZBWLTLSemvqRTdbh8t8EnCqA8oqL9NW73W2Z2o30xuJh28e32K7lYrBvK3JomP5hUx7mqYv6cZA1du3EkiTYzoCS5Puz6J3tes52r8S6uMLXa0y2bjr4xwJtsnzV00aJAT0lVeOW0RseB8InHHy0HanHCfRU48/mn0X/46jpwHYCuuHHA0YU3RaF4oqOiasGbPqNoaTUE1BPUT0vQVEU+/wAPi8oTTjXHhEcREI0T5OdqKvah/NB5Xj5croGmmmgaaaaBpppoGmmmgaaciP4jnqem2ikYtAhmScdoCnJCg8uECKXK8IvyX5agPINj72jjMv3N52ZCbAXTriZdYaV8zUg9YCN3hpvhv1BVT5Hs7VCOmuKC4hcuCbICbYsm+CNxpyE2ZeZAkZRjPtRxsBZeRQEkU0VPiiyDr4EjSq9IFxlQjMihG+LjbjhiKkocRmhBV98RFIlRtFjijvwBHTXUZ+FWMQvtOp3jJcX0W4oqY+DQUd9b4FLkuR4VE8eedQGQ04apHdjPq1IQEaF8kOXH4NFfbVWvgRD7F7eC8F8/qHdp9f3/AD/f/X/BNS4LSErEuYcqP7lCclNyXmZDDjgSWJzURYSCTgM+o244rREb4OKoqpNCf4erV7g7+7UbZ0m6uRZfmlRUVuzVVX3m4TRnIesMdi3kRZ+PwH44x0B6dcw+XYaNu9qqCi6TfdygY7+0usYNN0Tb13E6TEgtVIbfWXv8gFVyE3C3UwaRJlRTFPUadahNyFcNtUP3X1wRVQlEtjWx29VhvhEyjNImEu4ztNLfqJG0uWXT7jFnubXtxHmb3OI9P7oo4/i0ue5GdwqQM6xdyqisWreWxRvAsAtMNHurA9pBs5lG9lA/eUXSZRYzl83EsUuoi02T7qZ7itPYunZZ9AF2TFr8Mw+/hq/jtAEi6Zy2waosrcm0cimarpO2boJyfJs26KulLNMpurO6yzKenzam8ynJ7KUcm0tLaZhVK9PtLmaaeo5YWMp9x44i94ySdckHKEmRB0MsuE+ifNV+X5r81/rX81+eo6aaBpppoGmmmgaaaaBpppoGmmmgacJ9NNNA/Pn8/r+enCfTTTQNNNNBFPKpz9U/9+tf/tMP/NPsP/b70np/YnVTs9wn9n5a2AD80/rT/wB+tf8A7TD/AM0+w/8Ab70of/up2e0H0v6aaaDGPqp3vwzpo6ed5eoHPJdzDwvajDb3LcgPGYbc++ajwoyg5KrojkqC3KtBecYJgSlR0a7TVHVVeNaBOkzCM5wDYzDU3Slw73ePKSbzzeW+GW5annW52S8vZLndzdvssSLbIZrTFd79YSIov2pACOkx7oCubCfbPXsxelLDdpHXmWsH6pepTaLpi3Ujkqty02o3IHJ1ydmhmp3fYOQurUQfdrkGpSxu0vwHO9OMe4MUYMKJBbFBZhRmIjA+mjSoxHbFppDBCJO/sFO5efiXz450HcIGBkRvOSCJBb9V3hCRtrv9NOEUvK95d3nxwnlddmmmgaaaaBpppoGmmmgaePz5VPzRPmqfnx+/6aaconlfknlfPHhP3+eP6+F40Flzz7J6rdyXh95t1LrduplM7bYPutBsBlVBT6t1pnJse3EcfjQkxuybcn1RYWyy7cNXEccgceerFgNjMvKhfMS7lUXQYAjd9UpHqISx3hPhEUJYAbjXlU7RXkvzWls7wfFdyMTyDDsyqmb/ABzJqubWW1TPBLSvnU08mTZYsTdRn3z0DYE4hG0HubiEYo4pqiYYdPlr1c4fc7gbMZxtzByHD9tdw8fo9rNzr/PZ86wzHpzu4+RGlnbzncZVzK9ztuRqMVq8hr33ITN7KyR+cVpAWALEsM+VMEddYIxCQw76LscuUdBxOUUVREUfCpwvBL5VOOU86KaIbbRIouvK4DLSone5IZVEdiiiKqe8N8/GCqgp2lwa8ebWWF3u/Dw3OZNbhuLnldPf3MTb6ol5/Og02R4hDtGI2NZHkd0GLSnsYsLindkTXauNWXjdLIH7Nanz25RS2p7cUuav5VDnUuVVNJiTWM5ZQ3eNWuLtis7Kp8qmcwnIUyELQnKdulhQr5uTRt1k0L1Z4OOWED7PEZAVsRCIE4i94DGOURN/EiNATYEi/Je9CdBO1E/NfPjUVRxeexh11REzUWxFSQG1FHC4Uk8ApD3ceU5Txq2NZheUyMZ22gZpnttZZZhcumsb3MtvGywSNmGVVVfKhWKWuKsS7Rt/Drd2XImOYqdu4z6zEFxZprCDvqXK8Vby+gl4q5c3uPsWTIRAyDHrh7Fs4iC061I9WnvYzUx4pc5Y4jd2Stsm38DCMPjLI2QqFZsEZMeGU6IMmWy5IjME+2Lr8dtEI5DbZKjhMCJCpOoPYHcKGoqYotN3OdYljzeHvXV5BgNZ/fVuMYW4456reR39tW2dtCq61yOjzbshyvp7GQqmYMoEcuHVVR58tlgOC5DktRn+RYvQuZhV4rkWKw7+XDD16vCMxfqZOQ43DtO0nmqmwkU1MsyILCBLchRCNA9BNezD8IxPbbGaHAMMx2Hi+KYRWVtLjONVUZuqYgY7AjrGgx4EVn1m4zUVoI7RutkquiScthxwoSK73e23xz77rd5ZW1w7asY9Jz031f7cVYyzn7uOWqgySgNv2l7r6CPqvC+ogcanUjLhj5tT4V93Mpf+2MduchDLI9W05hkD7HnVUJKSzuVmC7GyC1+1Pe6evGC63MhV1m8clgowtvVITMP8UnWWfWMlMiVsZDER1xeXknySQCjskSqkLhl/wi+B488iMm2kYbEUAW09Fp+QRxDbbMAbdlctfFYekSlHJBXtbV9O7z5C3EXcp2W9Tt/yb7qD9pbm2u2r8dccisyIC1YWy/fSyRLckjbfWS1arW34k6/JCbBMqxpXiFud0srL3cmy2NbYtTVmI1MuhPCr6tyF6xu8jdchThyiJkMA6mI3ArKScsSNRvtzZ/2rFdOQ6xAJtGTrWQAEMoSJ4Y6i62DkqYQvTGAebBsVNGjX3HyPuwefw/HCcefOrJM8KrIRzbcGukFVn6kSDEUSJsi5FpVNxWAFVQU+a88/JQpSRV5h99Km6hZJXxsNg0tzEscNWhA51pdSZ1YdHaxso+0BcqYlZXMWkWZVhUyxs5E6NLKVGWALUinxx3dRXKw39x6WSkTcy6vpaycFZffkbW2K2hV+3EV8r9Fg2teT9Qs7NBB525KrccOjhrLUWbpaaCiY2EtR7DM565Flbh5jAp4IxBunWqjERq6h+rNzBatAIMZkzyeGwsZMd14plgwzJIRIE4l9Dt2lHNw+aWb7j3f3Pwc8JSFkuXSbmsybvdq3VyrKoTrDQ2+XD9l+kFw6aOi1YWIoP/ElxcbTQWvc2V2lPEb/AAItvMaPCsnv7LL73F3oDLtVY5fcXzeUWuQSGVAQKZOyEVtnSUFL3lUJTVfiW4U2qgWkGwrLJiNJiX8GVW3inDbNp2DNbUJ0QYZGoux7BFJuU2booTREPnnXu00FDxdvcQgbdLtJT0cHGNvixJzAo+L4+w3EpanFGKtaaDW45TNi3HrvTrkBAEHFGOTQGKuK2iFRPsvcnsbfpkcwaWSLX7Gbt7z9N+3vxKUiHtxshnlhgOABaDx2XFymLUMQLeebkVbGb6s9W2i/BW9ZtqYkqKiqAqXpm6rLLyL8CtPkgny2SFz29pcmg/L5pjv7PBsNu8t6v+muObk+JtdvUe8zNw5GGINpJ6pHsi3Xl0kqmF2Q3CiYU7LKjq7YZ8k72MI2DkCpM/dBDZnppyi/Jefy/tXnhP7eF/gv00RUX5Ki/wBXnQNNNNA0000DTTTQNNNNA0000DTTTQNNNNA0000DlBUVXwPcAqX5CpkgD3L+SKRIKL9VTnxrAL2lwGXSjMbUVB53f3pV7WDTh5Fa6q9oeUIfKIrjQLJbXuVCj8Hyiqg6z7Lu5RQHlee1TE1F1sD+Bwmw44e5bIhdaIm0JlTRTTnhcHOs6EG5+T9KHSTJkrDHqs6ksQx2RnSNIs/C4mw1fK6noSVWM+ojORxr9/ZuNicsnbqnWmC6fsWxsCrghyw+kXTTTQahvbJ4k1O9nxvJumzBlv5x0zR6rqO2nZjNuOsjuRt0k77trkVG0Jrb0CJbzFnUAyAGV+H/AMU2jXK4s4/YLa0FJYk8j5zKmBIceAwcacdejNuOEy42ZgbPeSoBIXyRU45Tzv0mwY1jHnVc/wBFYtvXFEVpARt6QyTT7MsFFUMQEQkB6fKl2KSr518ytPt+XRDvNN6M7FqTD2OocNqcn6Ydw8qdSumZfiZybEcu25WS568a+u9nO7G1yHKZE+DIvEzasFulhpEJXAyN01FG3HAdaVt6O5HWQDb/AGA85LThtWSkxVMAgSWk7vVitvShJHBVHvh88XDFRM2Ad9JyW1H7pgpFcgt8F6r7jQq964L8PHxh28L5886COmn9vP7/AK/v/t00DTTTQNNNNA0000DhPp+7+z6aaaaCHCfRPp8vy+mo8J9E/wD8+X8NNNA4T6JpppoIcIvhUTj5ccfl9NR4T6fu/s+mmmghwn0Tz8/Hz04T6J/DUdNA4T6JqHCfRPPz8fPUdNA0000DTTUURFVEUkFFXhSX5Cn/ANZf3J810ENNAIeRV8HhBHpMcfcwSWUpAI0YlohKx6LJCKKqIpp8Y/F9Zfa2dbj9VJvLywStpqmG7Nuree23Eqq9pltSMp0z1nFitEiKaOCy8XeIMi2ROIqB7SFXjbitiBypZK1EEx7kJ4BV40RF8c+g08qKqogrwflRQVwI2A6gdntrfaLdTe1WY5jEYyzqTpdnLjbmJAZdsGRl7V4ceL3eH5LMIY0CLllrItftXGY5yzS2oaq3myDrpLLcB7J7ZLC8r9oZZOQdqczyHBek+r9aNuJv7g9jPosw3Gtm3gBnANlLMo0WTVMwZQPP5buCZ+800umbxtjG7OPkLtrW7lcb6JOm/GOny16XIeBVdrtTcV8r74MXbbFtk2U5RZWMa7l51k105Hact84s7+OmVWGSONjJk5KAWpNo6nZoMbDYPvUF9ZAbNjlllsYbhTH4jzjDbzIG8JJDD1GXk717S8ePlrgyvLQIoMAYIjboxh7W/WBEF1VT83FcQlcX8zVV1jnuBtFu70Fz47dFX779V3TJbFImfeaZOc3Q3+2ZyMm3LPJbfNbGykVi5ftnaMs3F4/lBy6iZgoxKzHItDft2rllCvLiOa4tn1BTZXhN7ByPGMgh19rAu4TrJMjU3NWtvTWdi22469DK4iFHkw2HgR5yJIKS6LTjax1CrNNdYERtxnC9JCdBz1AYdWQyJtn2cg8TbRGJ+SBVaDuHzwnyXs0DTTTQNNNNA0000DTTTQNNNNA0000DTTTQQVCJUATMe9F7vTTvMhBO9URtVHvFe3g/jHsDuPglFBXD3I6GfvH7SPoQwzDCYbuenK3z3qrz1LN824A7ZWu2WfdPrI4+YtOLYZOWc59TK9UOtwo7dGs+0SxN6MEF/L8z7HGCJz0gJ1GSNFVXO59PRZAGeER5HHnGwcQnG0bbI3kUlbRsrP8As+63+Ujro65t18kjuu3uxLe2PTDtuUQvQoh22y7CcK3py5JyAjjNreBuR7z6Fu24BwIKLSOR+8SdEN4mmmmglpRX3fJkjSOh+M2C9ygfz4ad+FUT5pz2fLjx9LJb89MmyHVJgh7bdQm2GGbr4Q5awsgDGM4po2QVkPIKsJAVFxEblCgtTqn3uScB4RQmTeMk86v/AKaD5192ujLrL6Vzds9rlb60On2tk2E97D7V+Rj3UntxhdarChTYpNRvJWuo/O8gZf8ATZkXMva0a8qRoROUlqawrEYR1X7H5xmcva9/Jn8V3opm7GJkezu5LDuPbk4S3WlHGwq8zq3/AHmpr7qqcfabmxY15OEScH0nnfKp9RlsPdE4JE9L1QV5SX8MGvi7idDhfVbReO5pVFC5T4k41YPe7pt2K6nsHd2z6hNo8B3hwwrSuyJ/D88oYdzRyLmuGSFXayYclowCeyMmWkZ7l0mUcd4Uu9dBp5gWNfZRFnQpsSTCBw2Slx5LMmMJt9vKK9HceDz3fCnKl4XuFPHPsVeFISRQIV7SEvBCSfNFT6pz+Sqn79X5ufYvdJJ2qu7PXnUP0rYuz2NLtd0v7wTtl9shsG+9LDLYuIU9POiPXt4hRws7UnW3bAIMRHBBI482svPZNdRGPzRr+nnr6yjEdvEjtPBUdQezzfU/n6XbncllJHdG23KwGW9QvCMVKuiXHwbqiCUYSn1lkjYUwpeOAEnTXwDTfBOukvyBsVVO418qicp4RV58a7SAhJRXjlU5DheUdT/7pU5Q/n9U12P+zC9oDXNlZRvaIbbW0qvRZcWvTogqYTkp9oS9OO1MTfmQUQ5CqrCyhYeKODpvC06oemVj02O9sTF7mWej3o7RiGKjFjxesnLUJseeG32Yf6OAgBqKKpMo7/2uPU/DRSC9XC/u/in+euBEorx2kv7x4VP486s9/IZ7Zr9kDo2/t6z8u5/t/wDJu1TdpkHWZtdKXFd8ugTfTIs6BsJ7tj0n1Uff7aVaqZytcyxuDfObWT3chaRt77ap1xQGawijC1YTPWJWwyFQiXwLbir+SIicr+5Pi4/iqJ+/UFdbFOTMQ5ROO5f1jXnuaTjnlwOEQ0+SKqcEv5Y2DvJvYS9n83R7Q3tNO0xDYHHhEx+ag6n8pf4jKqid7a+C4Tn5aom/66NhNtJYY11AyMx6XdxDT3iTtDvphMnF8yZx1/j7GvJlLUTMhhfYVmAyPsQksVcUWZXc0347gzINxG2ydITRG22XHR4RTZSQhK0LqIvAmSAS8cr4FV512knaqpyJInyIV5Ek/IhX8xX8l4T+rWIGJ9eXR7m2S0uF4pv7gj+Q5XKiQqav7bWExYm8Dy1yy50qtZjxREUeQzdLtFXB+usm4GUYxMKPGiX2POSn+0G4kHIINibrypyox0E23T7uFUWxZ5FE8c/PQVHpp+ap+aLwqfmi+fCp+S+F/guoIqKvCKirwq8IvK8J4Vf6kVURfpzoI6acpxzz45Uefy7h8KP9aL80+afnpynKpz5TnlPzTheF5T8uF8L+/QNNQVUT5qiccKvKp4Rfkv8Ab+X11H8+PzVeET8+fp/XoGmmmgaafVfyTyq/RPqun+fH9v0/r8L/AA0DTTlE+ap8+P7fp/XqCKi/JUX+pdAVePkikS/qgPkzVPKoCKqcqg8mqcpwIkv5ca4qScKiipcgrqAiIquM9wgLopzwTTneCgvKKomiqieePBbXNVj9dNu7uXFr6qrivzZ9pJlJDSsissm45MbfUD7DEBUSFERXWSdbFVU0RcfML3C3t6pWokXoa2lf3Wx2WZi51A57Lmbd9P8AEgx3wq8l+5OatVOTWuV7gYxZSmkf2+l4xQQJ7US3D71R0hgUkK63V3v2t2XgV9xn+UxMbYyCwaxiBGjuvuWd3bFHkHCxengx2DdevLKRESMw2StxPWQgcmCvZ33m6cehzfPqVu6zPesXAKDavYOqdhZNgPTjXZBIvch3SbnOtXOLz+pFmTTU8Whk4Y03Hbt9pWWsqrpOXFAvEy1lzFmI9lnT0sez7276fraTuNk+X5nvjvxZ443QXO6e4cj3lMfgyX4NlkuP7QUD7tim12AZRe1sC6cwyBbW8doqypaWe8sAHD2GxzVxvuUkLky4Ts9NQTnlGyHlfiD9Ul58qnPjQStmuVlhhGYzEUmGxQWYxi2C96J6rRGDQorYLyqL6aeoYiaoipxr0pFfb5UDFz0wQWQNe3lfCKTp+e9URFXlRTuX6c6mWmgk4VziC93k0XqoaEHpJyXqr3PCrikv4bp+VbVFQR5Dkuedat93vZe7WV9nlG7fSdOtemXeyxtrfPLGvwCxKk2a3i3Hspj9lMut99v4IR42eWE+PLtqavvZk5p/G3LcblmLYHAGDJ2wa8T3Ku/CqiYoip2/F3D2+fUD4U7UXjjlV58LxoPn6xzf/c7BstrtsesbZmXsHuVfzGYcLM6O9dz3pxzC8k1Um8x/EcH3jsanD59rmrmHxLC6vKSfhFRDr3ai4YbtZhRGSlZVtSoz7UZ9iQy+xNYZlRH2TR1mTHkNI+w8y6HcDjbzRCbZiXBiQqKqiprYFuxs9t1vxg2T7W7x4VQ7jbaZXDj12S4blsRq6xnI4bcyLao1aU0hpWZCxLOFFOKJOfhk02aKnb2rqszz2eO/+xLEm+6Gd359zisQ3JK9JvUHcS7zA7GPHeSHjmH7X7mIj8jp5wbEKV30qmkqMCzZl1iqrqxCjNOe8MhdYi7CVshIXfRCQDKonqPMm4213tIi8KguutgXcQqhEiIi/PXNU7SMF8G2ZNuD+YOAvabZfQwJFEkTlEVFRFXWDznW5t/tXaQ8B6t8fyfpS3PtJzECJjO6NXIhY3uRew5jdPaW+0ORwRmLuNgjeROhFh5LcVmHyJsWTCsnKeMTpxm804FhDs2FmwrGFZtOuu+rJrzYchNykNRlRWHY777TrcV/vYF1DEnEBCJpoiUED26ahynnynw8IXn5KvyRfpz+XPz1FfHPPjj58+OPPHn6efH9egaacp9U1DlE+aon9ugjpqHKKvCKiqvyTnyuo8ovHC88qqJ+9U+aJ+9OF5T5pxoGuDrgsoJPIrTZ9vY4acA4iuC0pCqcqogZcGqonHC/Pxzz5TnjlOfPj8/HPPj93C8/TheflroJxuOkmW656IBHQHVVRIXoKG27Kc5MmwZFoQLuNSUUQVJVTzwHeqGi8ek55JRAkQe0+FXggXu8iQIrwqvCq0ndxz8OuKGCuOAJioMGIyX0XhiMht94OPGqIqNkqg0hAJr6pgPHCqSYiZt1q9PeL5vabOUVxP3h33h+5rV7A7U0i5XuXmzc+KxkENcahFJg0925R406t5ZxnriGMOJVTZAuOPRRjuTzHNl/ak9QDjcecGy/RbhbQv3uKZhZSZnVFmmRVc90iqsXyXZ+3q9oa3a3KItbLZl5I9AzrMErLWBOx5pJzL32qIZMy5TEBh2VNdbiQmUT1J0l5iPDFw0QmGSkPutti7JEgWOJKIueoCdyESJqQrmmIjwjmUY40faKk25e1CONqQoSg4KTVQTDntMeV4JFTlfnrjR+yC2QuVjTd9t4OpXqSgTWwe3B2j3P3Usch6dc9u3CFyc7J2XmxHquFjELIO3IMZxZbeazicqBUBHnz1qGnXq4Y9i17KB8Fdd9nx0uK4rjqEre2NK8JKLhCpK4rIKpFxyaKPwkqpyvHKhj3ubvtjuFQqaqxQV3I3NzuzPF9r9sMRkRLa7zzLThSJwVypFfdGtpqmvizMjya7eF37FxuptbSPDsZURmvk5sez66VbfpG6Zsa2lynIIl/lM3N90N2MukQI/FfTZNvZuXku7V1hrEsnVdu6rCrXL5GLVF6+zAdvIlREu3KmoOUtZGqnYX2fvRZ0n3mQZh009Lu0GyuVZHRLjl7ebd4RW0drbUgzY9kFVLcjIByoY2MOLL92JwBV5hs1L4URcrYLDgPgQqgIXu/rtdnuryuDDEU94AVcE+0URG46kgtN9qI4qgiKFTaaaaBpppoGocJ9E/gmo6aCHCfRP4JqOmmgahwnz4Tn68ajpoGocJ9E/hqOmghwn0T+GvJIgQJBE7IhxXnO3hXHY7LrnaKLwiGYEXCeeE54T8tezRfKKn1+nz0Fmty9otrd5MByrazdDAcfz/AG+zankU+UYpe1EK1pL6ldICKDbRXwbB0ZnAkTRCQqjRcl8KJrX3aexc9mosB97C+lDabZvMhFCg7sbK4/A283UwawL4vvFgmaV8OVNoL+IgnHjzWGCNWJUhOQ7uF2vpDYT1F4JTdEQdcUuXHBHnsQzVOS7OVQeV+FCVE+eu4WhAQAVJBBEREQvyT5c/XQaTZvsdcQiQZN3tV1kda1BnbLgyMYv9zN3LDefCWrAENGZmV7c2Z4vDzOuBkngCrl3dUAE6heuqpwuPmd9F3tQNqiq29m8r6aesKNepMK8sN0St+j+xwRYBsBEpaetwmh3/AGMyDJhkPzLedNnY8VG/Rw2I0a2Gzdeg/RkMZkEcFsEbF03HHRb4FHHHV5ccJETyZr+sS+V5X666Er4yCgCJC2Hoqw0hKjcZWBMQKMHyZJRMkJR/WTj6aD5fLXIus/aQn8a3+6Ed5LvLofqXM286TCrd6NiI1TYqkn313cDMZe0uRuZHEjMuysmEcCUKxxr0okiy9dTCSYJ1n9L+6eUP4ztxvDjeU3rMGxswgvfbuPE7RsyIjcmznXmQU1bFflA4/ESRVIpG4493tvGLJlr6pCisE4DqgnqghCJp4JANUUm+75q2SiKkPyVRHn5Jqym/3TL0+9VOGMbc9SOzu3292BRbuDkkfD9ysZrMsx5m+rY82JAtm6y2YkRhnxY1lPZZkoCOtty3xEkRwkUNKdfPpr5pqRSz4mQRFckCM+usY8tt5XTFxI0hI5PRnI0H01CISSVcjNmYKynqePechCedcSVHkH67MtHGXHWjiI0DrTgtKLBi4MgnwIVRURRbXzrI7LvYsdHjb6Js/ab79L+LEjTMHabpe3oyDYzaCHYOiS2VxBwDEYYU0e8viAX7yxBsXrV9lp2SXc2OrQTPZV9VVU/9lbYe0jtcS23pbBY2D4pnnTNj26+YUWLV4Pxaeqyvca53aprTObSFEdbZsbqfVQXrOSPvjrDZr2aCmNNSXMugP2kG1tQ3kO1vUXsp1g5O9MSplbXbn7cNdKWNVlaYOuS8vj7l4VI3st7G1rJMaLXx8Tdw+LBsWbaTOcvIrlWzGm2qDYj20KEyi9LfQ2ZKiigH1i7lstPqXnvJ8emFwmHGlREBkG3QfEjInW1bESC9Tyr2K2ip+N+GoE6TIOIio4oOGDbi9ioHKoo8KqJ+fGoCQyAkKDrUgVkHNFiP3etCHkkBG1cBvuGMjnpuEKJ2qScJ+S+ml9mr165ZRxr6+6/MI2hya8hR5d5tvifS1ju5NDgdvKVuTZYhj25tvuNhlxmldj0kTrY2VWOIY5Mu47IzH6avOQUduvKz2ReVZa69M6muuTf3cu1q3GgxO22EfndIkGkpbDl29hX9XhOWZq1lsywnM1T9bOlvwzrGYkpltt1Jhq2GKe5HUFs1s9VHke5O4uN41VNz28fU5MuRbk7evNvupCcg0sWxtYM2U3ElPnYu1/az6Dkd4QOSPElxzKeqzf1KqB0pdKObnDu2ol1Ub1dQrTG1PT7fbdPiLH3swHJ6E9wMoye6lsz6+6xKhusJxdq7qm7ArCzpHmwZd3adOvs6OjXpZyodydn9isIpt5pWOyMcyjfidS11jvZn7Fi/Bm31nuBuMcRq+yu8ye0rYdxktpYvE/cWzIzpSE8iKmaBwmHHBdJCV0VVRd7vxERfmCFxygL4+FPHhPpoNNWBeyIw6XZ0uU9TXUPvF1G5PXzolhb4qEn+TbpxzKJUuJIxqsybp3jWuWUcoqhxiFON48pc+07qvi3LgRHWhY1t5hKDIo0y0DLQpFRWxYCE3HF1lTNWRYV0VedcQVcaTtEF7hR0+OVngwY4oaKJOCaOCQul3j2OEhK3wqf9WPCIA/IRThPGjMJiOKAwhMsj2djDZdjIIAqKCDaIiCK88kieFJEX8k0Hq4T6J/DUeET5JxppoGmmmgahwn0T+Go6aCHCfRP4JqBIPBconyXleOfHHn+vx+WuWmgozJcbosop7WlyKqiW9RcVc2osYU1pp2FOqLGG7DsYkoHGnF9B+E88242iEJoXC61USvZGbV4k27adJu7W63SOtdHA8H2y2/vnLzpcp57vYM3IbXp8V7GqbJJuStOS7a8YfyCuR3KZo36yHXogsPbhPQbVF57lVeeSUuSVF/JSXzwieET5ceNcPdGUICRCHs+QivaC8Iop3CicLwKqifRPGg0C5ltF7V/aS4k49iuDdL/VdglS1Ds7Ld28z/JOn7cjLogNN2NpUwdlqDbvdHHmr7HCB2oxkz3MJvMHYsadPdxw5zkeNQlx1z4rgTZ23Ursj1N9JmMWLaRaLO+ozb+txzEszy1Xhlu4dUyMHyrcGwO3WoYtL9tZVZGiuRKiSyUgHTBF+jRIbHb2mKuJ3Ev4i9/wkaOen5T/AKsSQVAPkPanHy1z93a7kJEVOFVeEVURefKoqfmKl8ap+Zoi/loPn72u62+lHeazsKbbTevE8ktKasXJbSI5HtMTajUazI1ULzM7IINfHkSUmWMRr0mSKQrZGqNekJuDfM9yNuhXlc8xAkMRcBXMsojL03BQ20UlsEVUQFRE8J4T5a2K9QnSd00dWWOU+IdTGxe2O+uL49djklHQboYlU5fVVV+MCZVDcQYdvHksx7BK2wmwUlNiLqRZT7Pd2OEi4ifzKfskE+Xs5ekFP/yQwn/+L0GPWa9QWyG22LXOcZluZilNj+PR25U63YtoVuteD8qPBF0oVS/Os3PVOUMZtIcKS4jz7amAMo483jY17TbobmPw4WOb51GV2k4Q+wMfxelzJcmt7yYSNQ6Krq5uNwoEjJbac8EKO3JtocU3Hz9aayKqWtrO1fsq/ZvbHZ9jm6mzfRH02bYblYhKfm4tnWEbU4rj+UY/LlQpNbJkVNxAr2pcJ2RXTZkF42XBJyLKfZJexwkXPNIUdHTeQV7y+SqvKNovk0aT/sI4XxOIn65KpL5XQfPEGXe0U3fCTadNfQ1Fw+rpBKJkH6ce403pxy2TaSJASULBabbrEd+ot9QvQPXat72fbUcyPPcKA1UyWHlljkXjHsyM73Ygwsh6uOqDda4s7KNDkP7QbB3kraLaKkobZGbTJdpsvCNKuD3sx+HNcexUc2uanDpuTYwy6MvGKpbV6NE3HFXRid9dENt5fButF2G6CB2C26SJybYJ2kAL4EwAk8gnHo9Bte3u7j7WvS+MlLuFeEVS5/WIkTgiXyvK8/PQWU2F2K2Y6e8BqNr9jNtsI2x26xx+3docTweljU1JUP2trLsrRIUFhoQYKVYTZciUYH+NIedcVPxFTV8OE+ifw10R4rUb1EZ7hBwkJGu78JtUThUaD5NoS/EaJ+saqS+V16NBDhE+SIn9mo8InyTjTTQNQ4T6J/DUdNA0000DTTTQNNNNA0000DTTTQNNNNA0000DTTTQNNNNA0000HS+2roiKKiJ3ip/kvanPPav5Fzxwv0514ThuuOBygD6bjxhIIlddFCIextBJBTscTkj+NEAgDhC+aTTTQS9YjhEKq68o8EBITxEhNuKilynj4k7UQV89oKYp+svPnkVYSQab7QijHMUaJhE9QWQFURtlxFFYxKqB8YKSoKKHHC8pONNBL3ohvNdneqKg9qoa+ojop57TVVTgiVE5PhVRFVERUVdcAiu9nkAZUQY9JltxTYAhHtMUFQFOwfkK8efnwi8JqZ6aBpppoGmmmgaaaaBpppoGmmmgaaaaBpppoGmmmgaaaaBpppoGmmmgaaaaBpppoGmmmgaaaaD/9k=)![ref3]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEADEDASIAAhEBAxEB/8QAGQABAAIDAAAAAAAAAAAAAAAAAAIFBgkK/8QALRAAAAQEAwcDBQAAAAAAAAAAAQIDBAAFBgkIERIHExQWGVmYFSHWYWNxktL/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A639oOAQdpNbVJW441bgFBjV00ezs9J7NsR3K1E06dXScZbSsh5NfejygmrSmy4txpKUob0csxwiVWyQWlzNUbgtzsgnRKbQnivEiZc8/YheQhyL9Mx/MIQFh0xC9wi595Yj8Bh0xC9wi595Yj8BhCAdMQvcIufeWI/AYi3tqCwm8uXC4BczeEaqFf8I/xV8SzcnZrJKFQdoDQZN81WABScoiYu9ROcmoueYIQGxr1Nz9v9B/qEIQH//Z)![ref10]![ref10]![ref10]![ref10]![ref10]![ref10]![ref15]![ref10]![ref15]![ref10]![ref15]![ref10]![ref10]![ref10]![ref10]![ref15]![ref15]![ref10]![ref2]

<a name="br20"></a> 

20

Fast Approximations and Coresets for (k, ℓ)-Median under Dynamic Time Warping

Algorithm 1 k-median framework

procedure k-Routine((X, ϕ), ε, A)

▷ A is (α, β)-approximate metric k-median

~~p~~

a ← Θ(ε

−1

~~log(~~ε <sup>1</sup>~~)~~), b ← Θ(a<sup>2</sup>)

▷ Determine the success probabilities.

−

√

s ← a kn ~~log~~ k

Choose a set S of s points sampled without replacement from X

C

<sup>′</sup> ← A

((S, ϕ| ))

S

Select the set M of points x with the b log largest values of min<sub>c′∈C′</sub> ϕ(x, c<sub>′</sub>)

M

kn

k

s

return C = C ∪ A((M, ϕ| ))

′

end procedure

procedure k-Median((X, ϕ), ε, A)

▷ A is (α, β)-approximate metric k-median

~~p~~

a ← Θ(ε

−1

~~log(~~ε <sup>1</sup>~~)~~), b ← Θ(a<sup>2</sup>)

▷ Determine the success probabilities.

−

√

s ← a kn ~~log~~ k

Choose a set S of s points sampled without replacement from X

C

<sup>′</sup> ←

k

-Routine((S, ϕ| ), ε, A)

S

Select the set M of points x with the b log largest values of min<sub>c′∈C′</sub> ϕ(x, c<sub>′</sub>)

M

kn

k

s

return C = C ∪ k-Routine((M, ϕ| ), ε, A)

′

end procedure

with A is a (3(1 + ε)(2 + α), 2β)-approximate algorithm for k-median in metric spaces with

constant success probability.

▶ Lemma 33. Let X be a set of n points, and let ϕ be a distance function that can be

computed in T time for any x, y ∈ X. Let T (n) be the running time of the (α, β)-

ϕ

A

approximate algorithm for k-median on n elements. Then k-Routine has a running time of

p

<sub>O(n</sub>2<sub>T</sub> T<sub>A</sub>

\+

(min(<sub>n, ε−1</sub> log( ) log(<sub>ε−</sub><sup>1</sup>)))).

kn

k

ϕ

Proof. The only steps that take time outside the two calls to A are sampling S and construcing

M. Computing the values min ϕ(x, c<sup>′</sup>

) for all x ∈ X can be done in a single execution of

c<sup>′</sup>∈C′

Dijkstra’s algorithm, starting with the points of C ⊂ X at distance 0, by adding a temporary

′

point with distance 0 to all points in C and starting Dijkstra’s algorithm on this temporary

′

point. This takes O(n<sup>2</sup>T ) time, which also dominates the time it takes to sample S as well

ϕ

as constructing M from these computed values.

◀

▶ Lemma 34. Let X be a set of n points, and let ϕ be a distance function on X, which can

be computed in time T , and further there is a constant ζ such that ϕ ≤ ζϕ. Let Y ⊂ X. Let

ϕ

ε > 0 and let A be the (10 + ε, 1)-approximation for metric k-median of Lemma [30.](#br19)[ ](#br19)Then

k-Routine returns a (3(1 + ε)ζ(12 + ε), 2)-approximation of k-median in the metric space

(Y, ϕ| ) in time O(|Y |<sup>2</sup>T + |Y |<sup>2</sup>k log(k)ε<sub>−</sub><sup>2</sup> log(ε<sub>−</sub><sup>1</sup>) + k<sup>7</sup>ε<sub>−</sub><sup>5</sup> log<sup>5</sup>(|Y |))).

Y

ϕ

Proof. The running time bound follows by Lemma [33](#br20)[ ](#br20)and Lemma [30,](#br19)[ ](#br19)together with the fact,

p

that min(|Y |, ε<sub>−</sub><sup>1</sup> k|Y | ~~log(~~k~~) log(~~ε <sup>1</sup>~~)~~)<sup>3</sup> ≤ |Y |<sup>2</sup>k log(k)ε <sup>2</sup> log(ε <sup>1</sup>). The approximation

−

−

−

guarantee follows by Theorem [32,](#br19)[ ](#br19)Lemma [31](#br19)[ ](#br19)and Lemma [30.](#br19)

◀

By combining the presented subroutines, we obtain our two main results of the section.

The ﬁrst is Theorem [36,](#br21)[ ](#br21)which provides a linear time approximation algorithm for k-median

in metric closures, assuming the underlying distance is reasonably well approximated by its

metric closure. The second is Corollary [38,](#br21)[ ](#br21)combining Theorem [36](#br21)[ ](#br21)with Lemma [26](#br17)[ ](#br17)to yield

an approximation algorithm for p-DTW with an unoptimized approximation guarantee.

![ref3]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEAhQDASIAAhEBAxEB/8QAGAABAQADAAAAAAAAAAAAAAAAAAYFCAr/xAArEAABAgUDBAICAgMAAAAAAAAAAQIJGVmY1gMEBQYHERIIFhMUFSEiYXH/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A63+4PwEXuT1t1J1uvzViAdBr1dym95t/Sfbb5HfVuienX6vq9eN6V4H6bvv4fiGe3rp7L9vcerWtT8q+PKxHFQyU1uO2eqsQWJ2xX6LXemn8r1Zpt8+f6Y36Evhv+vK/9AAyEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAk7u4bL2c3pcWvz+iYaug/jtzu11db5U++5R2lq7ViaTdb6InjQf+ZXamn6/wCT2aa+yevhQA3x7FdlU7I9BaXQid2e9fdv9fl+U5Jesu93Xf33rzcLyGqzUXZbrqH+K4n8vH7L09OP236bf1mPe33f7eUAAf/Z)![ref6]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEAC8DASIAAhEBAxEB/8QAGAABAAMBAAAAAAAAAAAAAAAAAAQGBwr/xAArEAAAAwYEBgIDAAAAAAAAAAABAgYAAwQFBwgJERIUExYZWZjWFUEXIjL/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A63qhWCjUlZqdbjeriAIMVbMY6eHSdNrjuVkSnjvhKYZalpDybHfESh3q0uoLdxGgpShxRyzalSrDIB7LoJ6OINidkF5DOjaHd1+h2XMgDpIUEF+pQ+gzHIPtjGCf0xC9wjE+8sR9BZ0xC9wjE+8sR9BYxgdMQvcIxPvLEfQW06jFlhqGVPTy0Ldve3V4juXKKBOjq5V6/ICFiN7Au3JY2JT/ACpKeLHwAiL+XRO7Ltn48TQf+WMYP//Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEABYDASIAAhEBAxEB/8QAFwABAAMAAAAAAAAAAAAAAAAAAAYICv/EACcQAAAEBAYCAwEAAAAAAAAAAAECBAUDBggSAAcJERUWExQXISMx/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/ANb2YVAo5kznM87jWrqASGM2uK58PKeW1R3VpJl48YSmFtlZh6au4hoh3WwkXtqLClKHlHbfELa9MgIrciijqDanZBiJoJ7IdWAlhkuIA2kL0IbSB/ChuOwfW+GGAslTtSOFPEyPsyBU3V1njzrGVk4GonOj5MltntXpV/LsTX1tl497H1vTMv8APG3QKFKfxfreVhhgP//Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEACYDASIAAhEBAxEB/8QAGAABAAMBAAAAAAAAAAAAAAAAAAUGBwr/xAAqEAAAAwYFAwUBAAAAAAAAAAABAgUAAwQGCBIHCREUFhMZWRUyQZjWUf/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDrexCoEHEmdJlnca1cwCQzTcpxy2eU8NqjuLSTLp3tpxTZVQeGx3o6QS613BbuItKUodUdNRpKVlkg+TYJ6OYLmdkF5DuzWO6rxI7JqHtIXgQ2lD4DUdP6xjBIdsQvkIzPvtiP4FnbEL5CMz77Yj+BYxg17AakZ7T/ADAszYWqWsTGo60lxMujL9QGNo4jywmEBSg44qyjpPGUbZL5ATwgSKO4e2p0XHQ3RHcXkMYwf//Z)![ref10]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEABoDASIAAhEBAxEB/8QAGAABAAMBAAAAAAAAAAAAAAAAAAQGCQr/xAAmEAAABQIFBQEBAAAAAAAAAAABAgMFBgQSAAgJERMHFBUWQSEi/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AOt7qFkEHqTM5POBzq6gEDGWuNc+HiXTbMd6tCI8dW0wtsWYPTa7xDQS+1Oi7uotKUoco7bjSmrTIBVuolR1BtTsgqUyRrE819iZdyANpChAv5KHwNx2D7hhgJNTpk8FNUrF1BtTwxkqdc5SnzYXEESonEAMUYF+hv8ANwxorBo+aPQmHsBn6RyAzHFo+zmfpM6eVkj2ZsaaSiF3kDpwIeSe3IUBrXWv4Ee8r1qio4k+SwrDAf/Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEADMDASIAAhEBAxEB/8QAGAABAQADAAAAAAAAAAAAAAAAAAUGCQr/xAAtEAAAAwUGBgEFAAAAAAAAAAABAgMABAUGCQcIERITFhQVGSFZmNYYJDFhcf/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDrftBuCDaTO0yTuN9WoBIYzdFH2NnlOza8dtaSZdOrlOMNlWA7NfuTwgmbKm5cW8ZSlKGqOGI4RCqZILQ5zVGoLU7IJ0SmyJ3rxImXHHsQuwhwL+sR/rGMFDpiF8hFT72xH4CzpiF8hFT72xH4CxjA6YhfIRU+9sR+AtHidNUHRRFD6+qlzyC6TwoZR7vUa6yYIaeKKSmxC5EV9T7hPAdXTS7lydzGDYzZXJw2fWdSdJIzXOc8jLMDdIRu+0SObmniYeFKJeaTRH+FcubRd4xzPT7wjvrGAB0i/hjGMH//2Q==)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEACYDASIAAhEBAxEB/8QAGAABAAMBAAAAAAAAAAAAAAAAAAUGBwr/xAAqEAAAAwYFAwUBAAAAAAAAAAABAgUAAwQGCBIHCREUFhMZWRUyQZjWUf/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDrexCoEHEmdJlnca1cwCQzTcpxy2eU8NqjuLSTLp3tpxTZVQeGx3o6QS613BbuItKUodUdNRpKVlkg+TYJ6OYLmdkF5DuzWO6rxI7JqHtIXgQ2lD4DUdP6xjBIdsQvkIzPvtiP4FnbEL5CMz77Yj+BYxg17AakZ7T/ADAszYWqWsTGo60lxMujL9QGNo4jywmEBSg44qyjpPGUbZL5ATwgSKO4e2p0XHQ3RHcXkMYwf//Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEAhQDASIAAhEBAxEB/8QAGAABAQADAAAAAAAAAAAAAAAAAAYFCAr/xAArEAABAgUDBAICAgMAAAAAAAAAAQIJGVmY1gMEBQYHERIIFhMUFSEiYXH/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A63+4PwEXuT1t1J1uvzViAdBr1dym95t/Sfbb5HfVuienX6vq9eN6V4H6bvv4fiGe3rp7L9vcerWtT8q+PKxHFQyU1uO2eqsQWJ2xX6LXemn8r1Zpt8+f6Y36Evhv+vK/9AAyEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAk7u4bL2c3pcWvz+iYaug/jtzu11db5U++5R2lq7ViaTdb6InjQf+ZXamn6/wCT2aa+yevhQA3x7FdlU7I9BaXQid2e9fdv9fl+U5Jesu93Xf33rzcLyGqzUXZbrqH+K4n8vH7L09OP236bf1mPe33f7eUAAf/Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEAGkDASIAAhEBAxEB/8QAGQABAAIDAAAAAAAAAAAAAAAAAAEFBgkK/8QALxAAAAMECAUDBQAAAAAAAAAAAAECAwQJEQUGCBIZWZjWBxMUFiEVJjFhZHGU4v/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDrf4g2BD4k12rJXc7asQCoZ1upR9ptdU+G1o7tapNXVtbqzo2qtA9mv3o9EIvXWbl1bxdSlJc05TPCKKhkk2o5zanEFidoNbFKrjO1eaGaZz8IT2Eck/SZ/kAAWGGInMIifasT2CGGInMIifasT2CAAGGInMIifasT2CGGInMIifasT2CAAGGInMIifasT2CGGInMIifasT2CAAGGInMIifasT2CGGInMIifasT2CAAGGInMIifasT2CIVDFSRGZRB4n0yIzKdrE5eCn59ggACnq9DhNbyzpJvb1iTPK6PpV4MnV4tR8xyfEOD2yVyX537GLqGL5cJL+i+jqE+Jo+RtL6H7x+/Y/gAAf/Z)![ref10]![ref10]![ref10]

<a name="br21"></a> 

Conradi, Kolbe, Psarros and Rohde

21

▶ Lemma 35 ([[29](#br25)]). Assume that C′, computed in the algorithm k-Median, is an (α, β)-

approximation of k-median of S in the metric space (S, ϕ| ). With constant probability

S

depending only on a and b, there is a subset of X of size at least n − b log k, whose cost

kn

s

under ϕ with C as medians is at most (1+ε)(2 +α)∆<sup>opt</sup>, where ∆<sup>opt</sup> is the cost of an optimal

′

k-median under ϕ of X.

▶ Theorem 36. Let X be a set of points and let ϕ be a distance function on X with ϕ ≤ ζϕ.

Let ε > 0 and let A be the (10 + ε, 1)-approximation for metric k-median of Theorem [29.](#br18)

Then k-Median returns a (11ζ<sup>2</sup>(1 + ε)<sup>2</sup>(12 + ε), 4)-approximation of k-median of X in the

metric space (X, ϕ) in time O(nk log(k)T + nk<sup>2</sup> log<sup>2</sup> k + k<sup>7</sup>ε <sup>5</sup> log<sup>5</sup>(n)).

−

ϕ

Proof. Let ∆<sup>opt</sup> be the cost of an optimal solution to the k-median problem on X under ϕ.

By Lemma [34](#br20)[ ](#br20)the set C is a (3(1 +ε)ζ(12 +ε), 2)-approximation of k median of S in (S, ϕ| )

′

S

and can be computed in time O(nk log(k)T +nk<sup>2</sup> log<sup>2</sup>(k)ε <sup>4</sup> log<sup>2</sup>(ε <sup>1</sup>)+k<sup>7</sup>ε <sup>5</sup> log<sup>5</sup>(n)). The

−

−

−

ϕ

set M can be computed in O(knT ) time. By Lemma [34](#br20)[ ](#br20)the set C is a (3(1 +ε)ζ(12 +ε), 2)-

′′

ϕ

approximation of k median of M in (M, ϕ| ) and can be computed in time O(nk log(k)T +

M

ϕ

nk

<sup>2</sup> log<sup>2</sup>( ) −<sup>4</sup> log<sup>2</sup>( −<sup>1</sup>) +

k ε

ε

7

−<sup>5</sup> log<sup>5</sup>( )).

k ε nk

Now observe that for the set O of Lemma [35,](#br20)

X

X

X

ϕ(x, C

<sub>′</sub>) ≤

ϕ(x, C<sub>′</sub>) ≤ ϕ(x, C<sub>′</sub>)

x∈X\M

x∈X\M

x∈O

X

≤

<sub>ζϕ(x, C′</sub>) ≤ ζ(1 + ε)(2 + 3(1 + ε)ζ(12 + ε))∆<sup>opt</sup>.

x∈O

Further observe, that there exists a clustering of M with cost 2∆<sup>opt</sup> by replacing every point

of an optimum with its closest point in M. Thus,

X

X

ϕ(x, C

<sub>′′</sub>) = ϕ|<sub>M</sub> (x, C<sub>′′</sub>) ≤ 6(1 + ε)ζ(12 + ε)∆<sup>opt</sup>.

x∈M

x∈M

Overall we get a (11ζ<sup>2</sup>(1 + ε)<sup>2</sup>(12 + ε), 4)-approximation (as ζ ≥ 1 and ε > 0) in time

O(nk log(k)T + nk

<sup>2</sup> log<sup>2</sup>( ) −<sup>4</sup> log<sup>2</sup>( −<sup>1</sup>) +

k ε

ε

k ε

7

−<sup>5</sup> log<sup>5</sup>( )).

n

◀

ϕ

We brieﬂy discuss simpliﬁcation schemes for curves under p-DTW (for more details refer

to Appendix [C).](#br27)[ ](#br27)We reduce the problem of ﬁnding an (1 + ε)-approximate simpliﬁcation to

ﬁnding a (1 + ε)-approximation of a center point for a set of ≤ m points, where the objective

is to minimize the sum of the individual distances to the center point raised to the pth

power. Note that for p = ∞, the problem is that of ﬁnding a minimum enclosing ball, and

for p = 2, the problem can be reduced to that of ﬁnding the center of gravity of the set of

discrete points, which can both be solved exactly. Furthermore, we show (Proposition [37)](#br21)

that for all dtw , there is a deterministic 2-approximation that is a crucial ingredient for our

p

approximation algorithms of (k, ℓ)-median under p-DTW.

▶ Proposition 37. For σ = (σ , . . . , σ ) ∈ X<sup>d</sup> and integer ℓ > 0, one can compute in

1

m

<sub>d</sub> such that

ℓ

m

O(m d

<sup>2</sup>( + + )) time a curve

ℓ

m

σ ∈ X

∗

inf dtw (σ , σ) ≤ dtw (σ , σ) ≤ 2 inf dtw (σ , σ).

∗

p

ℓ

p

p

ℓ

σ ∈X<sup>d</sup>

σ ∈X<sup>d</sup>

ℓ

ℓ

ℓ

ℓ

▶ Corollary 38. For any ε > 0 the procedure k-Median from Algorithm [1](#br20)[ ](#br20)can be used

/p

to compute a (72(1 + ε)<sup>2</sup>(12 + ε)(16mℓ<sup>3</sup>)<sup>1</sup> , 4)-approximation for (k, ℓ)-median for an

input set X of n curves of complexity m under dtw in time O(nm<sup>3</sup>d + nk log(k)ℓ<sup>2</sup>d +

p

k ε

<sup>2</sup> log<sup>2</sup>( ) −<sup>4</sup> log<sup>2</sup>( −<sup>1</sup>) +

ε

7

−<sup>5</sup> log<sup>5</sup>( )).

k ε

n

nk

![ref10]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEAA8DASIAAhEBAxEB/8QAFgABAQEAAAAAAAAAAAAAAAAAAAUK/8QAIxAAAAUDBAMBAAAAAAAAAAAAAgMEBQYBBwgAExQWCRESFf/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDWzcPAOtxZfKptXNfyBQWsscXOQDidusj+sQmPGKRVNq1Rdh6Yu/JZCdzbIQVVKNsoBYN0Xz71Rs7hCOzc8g87DmFnRdIKVM7EjhV48ge7wBfRyjK5HQxwjnUmvkmtwlFFzWPmA4jgnSqfRm18CaaD/9k=)![ref10]![ref10]![ref10]![ref10]![ref10]![ref10]![ref10]![ref10]![ref10]![ref10]![ref2]

<a name="br22"></a> 

22

Fast Approximations and Coresets for (k, ℓ)-Median under Dynamic Time Warping

Proof. Let X = {τ | τ ∈ X} be a set of 2-approximate optimal ℓ-simpliﬁcations of X under

∗

∗

dtw . By Proposition [37,](#br21)[ ](#br21)X can be computed in O(nm<sup>3</sup>d) time. We now apply Theorem [36](#br21)

∗

p

and Lemma [25](#br17)[ ](#br17)to obtain a (12(2ℓ)<sup>2</sup> (1 + ε)<sup>2</sup>(12 + ε), 4)-approximation of k-median of X

/p

∗

in (X , dtw | ) in time O(nk log(k)ℓ<sup>2</sup>d + nk<sup>2</sup> log<sup>2</sup>(k)ε <sup>4</sup> log<sup>2</sup>(ε <sup>1</sup>) + k<sup>7</sup>ε <sup>5</sup> log<sup>5</sup>(n)). By

∗

−

−

−

p

X∗

Lemma [26,](#br17)[ ](#br17)the computed set is a (6(4mℓ)<sup>1</sup><sub>/p</sub>12(2ℓ)<sup>2</sup><sub>/p</sub>(1 + ε)<sup>2</sup>(12 + ε), 4)-approximation for

(k, ℓ)-median for X under dtw<sub>p</sub>.

◀

6

Coreset Application

The theoretical derivations of the previous sections culminate in an approximation algorithm

(Theorem [40)](#br22)[ ](#br22)to (k, ℓ)-median that is particularly useful in the big data setting, where n ≫ m.

Our strategy is to ﬁrst compute an eﬃcient but not very accurate approximation(Corollary [38)](#br21)

of (k, ℓ)-median. Subsequently, we use the approximation to construct a coreset. The coreset

is then investigated using its metric closure, where by virtue of the size reduction we can

greatly reduce the running time of slower more accurate algorithms metric approximation

algorithms, yielding a better approximation.

▶ Theorem 39 ([[4](#br24), [21](#br25)]). Given a set X of n points in a metric space, one can compute a

(5 + ε)-approximate k-median clustering of X in O(ε <sup>1</sup>n<sup>2</sup>k<sup>3</sup> log n) time. If P is a weighted

−

point set, with total weight W, then the time required is in O(ε<sub>−</sub><sup>1</sup>n<sup>2</sup>k<sup>3</sup> log W).

Algorithm 2 ((32 + ε)(4mℓ)<sup>1</sup><sub>/p</sub>)-approximate (k, ℓ)-median

procedure (k, ℓ)-Median(X ⊂ X , p, ε)

d

m

ε

<sup>′</sup> ←

ε/

46

Compute (O((16mℓ<sup>3</sup>)<sup>1</sup><sub>/p</sub>), 4)-approximation C<sub>′</sub> (Corollary [38)](#br21)

Compute bound of sensitivity for each curve x ∈ X from C (Lemma [16)](#br10)

′

Compute sample size s ← O(ε <sup>2</sup>dℓk<sup>2</sup>(m<sup>2</sup>ℓ<sup>4</sup>)<sup>1</sup> log<sup>3</sup>(mℓ) log<sup>2</sup>(k) log(ε <sup>1</sup>) log(n))

−

/p

−

Sample and weigh ε -coreset S of X of size s (Theorem [20)](#br11)

′

Compute a 2-simpliﬁcation for every s ∈ S resulting in the set S (Proposition [37)](#br21)

∗

Compute metric closure values ϕ = dtw | (Lemma [27)](#br18)

p

S∗

Return (5 + ε , 1)-approximation of weighted k-median in (S , ϕ) (Theorem [39)](#br22)

′

∗

end procedure

▶ Theorem 40. Let 0 < ε ≤ 1. The algorithm (k, ℓ)-Median in Algorithm [2](#br22)[ ](#br22)is a ((32 +

ε)(4mℓ)

,

<sup>1</sup><sub>/p</sub> 1)-approximate algorithm of constant success probability for ( )-median on

k, ℓ

ꢂ

ꢃ

√

curves under dtw<sub>p</sub> with a running time of Oe n(m<sup>3</sup>d + k<sup>2</sup> + kℓ<sup>2</sup>d) + ε<sub>−</sub><sup>6</sup>d<sup>3</sup>ℓ<sup>3</sup>k<sup>7</sup> m<sup>6</sup>ℓ<sup>12</sup> ,

p

where Oe hides polylogarithmic factors in n, m, ℓ, k and ε <sup>1</sup>.

−

Proof. Computing a (O(m<sup>1</sup> ℓ<sup>3</sup> ), 4)-approximation via Corollary [38](#br21)[ ](#br21)takes time O(nm<sup>3</sup>d +

/p /p

nk log(k)ℓ nk

<sup>2</sup> + <sup>2</sup> log<sup>2</sup>( ) + <sup>7</sup> log<sup>5</sup>( )) and has constant success probability. From this we

k

k

n

can compute a ε<sub>′</sub>-coreset S by Theorem [20](#br11)[ ](#br11)of size

O(ε<sup>−2</sup><sub>dk2</sub>(m<sup>2</sup>ℓ<sup>4</sup>)<sup>1</sup><sub>/p</sub>

log<sup>3</sup>(mℓ) log<sup>2</sup>(k) log(ε<sub>−</sub><sup>1</sup>) log(n))

with constant success probability, and in time O(kn). Computing S takes O(nm<sup>3</sup>) time.

∗

Computing the metric closure of S takes O(|S|<sup>3</sup>) time and computing a (5+ε )-approximate

∗

′

solution to the k-median solution of S takes O(ε <sup>1</sup>|S|<sup>2</sup>k<sup>3</sup> log |S|) time by Theorem [39.](#br22)[ ](#br22)This

∗

−

![ref12]![ref5]![ref3]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEAhQDASIAAhEBAxEB/8QAGAABAQADAAAAAAAAAAAAAAAAAAYFCAr/xAArEAABAgUDBAICAgMAAAAAAAAAAQIJGVmY1gMEBQYHERIIFhMUFSEiYXH/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A63+4PwEXuT1t1J1uvzViAdBr1dym95t/Sfbb5HfVuienX6vq9eN6V4H6bvv4fiGe3rp7L9vcerWtT8q+PKxHFQyU1uO2eqsQWJ2xX6LXemn8r1Zpt8+f6Y36Evhv+vK/9AAyEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAk7u4bL2c3pcWvz+iYaug/jtzu11db5U++5R2lq7ViaTdb6InjQf+ZXamn6/wCT2aa+yevhQA3x7FdlU7I9BaXQid2e9fdv9fl+U5Jesu93Xf33rzcLyGqzUXZbrqH+K4n8vH7L09OP236bf1mPe33f7eUAAf/Z)![ref10]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEADIDASIAAhEBAxEB/8QAGAABAQADAAAAAAAAAAAAAAAAAAUGCQr/xAAuEAAABAMFCAAHAAAAAAAAAAABAgMFAAQGCAkREhMHFBUWGVmY1hghIyQyQVH/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A639oNgQdpNbVJW421bwCgxq50nXs9J7NrR3K1E06dXKcW2lWHk2e4O0EzZU5Le5jKUpQ1RwxHCGq7JBZtklRvBbzsgqS6Zsidq8SJkxD8SF5CHKUP0GI4f2EICh0xC9wi8+8sR9Bh0xC9wi8+8sR9BhCAdMQvcIvPvLEfQYjON2qEssWW+Pu8uXA0spMitM2qNVf6aqaYS4KciF+2MY4LGSw+ayaZ8wZcBQgNm9HsnLVJUvTgu75UAsFOsjIL/U7hxWpHvhTZKyHF6hdNGX4k+OWhvrs4aCO+T60xMaKepkKhCA//9k=)![ref10]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEAhQDASIAAhEBAxEB/8QAGAABAQADAAAAAAAAAAAAAAAAAAYFCAr/xAArEAABAgUDBAICAgMAAAAAAAAAAQIJGVmY1gMEBQYHERIIFhMUFSEiYXH/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A63+4PwEXuT1t1J1uvzViAdBr1dym95t/Sfbb5HfVuienX6vq9eN6V4H6bvv4fiGe3rp7L9vcerWtT8q+PKxHFQyU1uO2eqsQWJ2xX6LXemn8r1Zpt8+f6Y36Evhv+vK/9AAyEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAiWI2oRE+uxXAQAEsRtQiJ9diuAk7u4bL2c3pcWvz+iYaug/jtzu11db5U++5R2lq7ViaTdb6InjQf+ZXamn6/wCT2aa+yevhQA3x7FdlU7I9BaXQid2e9fdv9fl+U5Jesu93Xf33rzcLyGqzUXZbrqH+K4n8vH7L09OP236bf1mPe33f7eUAAf/Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEACUDASIAAhEBAxEB/8QAGAABAAMBAAAAAAAAAAAAAAAAAAUGCAr/xAArEAAAAwYDCAMBAAAAAAAAAAABAgQAAwUGCBEHCRITFBUWGVmY1hchYXH/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A638QaBBxJnaZJ3GtXMAkMZuii2NnlPDao7laSZdO90nGGyrAeTV3B4QTVpdot7UaSlKG1G1xpEKyyQfQ5G9HMFzOyCdyU2h3VeJHZb3+iF5CGxfy4/1jGCQ6Yhe4RmfeWI+gs6Yhe4RmfeWI+gsYwaXptpzPT6ee4QOPdSeOJI8olxWRfUPij8lxOA7ijiLoyaXFfAoLwxIvFSD2Iudm/wB6fJkp9RNjYxjGD//Z)

<a name="br23"></a> 

Conradi, Kolbe, Psarros and Rohde

23

is a (4mℓ)<sup>1</sup><sub>/p</sub>(32 + 13ε )(1 + ε )-approximation to (k, ℓ)-median of X by Lemma [26](#br17)[ ](#br17)and

′

′

Theorem [20.](#br11)[ ](#br11)Overall the approximation factor is (4mℓ)<sup>1</sup><sub>/p</sub>(32 + ε) and the running time is in

ꢅ

ꢅ ꢆꢆ

1

ꢀ

ꢁ

√

3

\+

log <sup>2</sup> + ( log )<sup>2</sup> +

ε<sup>−6</sup>d ℓ k m<sup>6</sup>ℓ<sup>12 log9(</sup>mℓ) log<sup>6</sup>( ) log<sup>5</sup>( ) log<sup>3</sup>

3

3

7

p

.

◀

O n m d

k

kℓ

k

k

k

n

ε

Combining the computed ε-coreset with the (k, ℓ)-median algorithm from [\[18](#br25), Theorem 35]

instead, we achieve a matching approximation guarantee and improve the dependency on n.

The improved approximation guarantee from Corollary [41](#br23)[ ](#br23)compared to Theorem [40](#br22)[ ](#br22)comes

at the cost of an exponential dependency in k, as is also present in their results.

▶ Corollary 41. Let 0 < ε ≤ 1 and 0 < δ ≤ 1. There is an ((8 + ε)(mℓ)<sup>1</sup><sub>/p</sub>, 1)-approximation

for (k, ℓ)-median with Θ(1 − δ) success probability and running time in

ꢂ

ꢂ

ꢃꢃ

ꢀ

ꢁ

√

+2

O n m d k<sup>2 +</sup> kℓ<sup>2) +</sup>

(

3

\+

k

<sup>7</sup> + 32

k<sup>2</sup>ε<sup>−1 log(1</sup>/δ

) <sup>k</sup>

md m ε<sup>−2</sup>dℓk<sup>2</sup> m<sup>2</sup>ℓ<sup>4</sup>

<sup>3</sup> +

,

e

p

where Oe hides polylogarithmic factors in n, m, ℓ, k and ε <sup>1</sup>.

−

Finally, combining Theorem [40](#br22)[ ](#br22)with Theorem [20](#br11)[ ](#br11)yields the following result.

▶ Corollary 42. The algorithm (k, ℓ)-Median in Algorithm [2](#br22)[ ](#br22)can be used to construct an ε-

ꢂ

ꢃ

√

coreset for (k, ℓ)-median in time Oe n(m<sup>3</sup>d + k<sup>2</sup> + kℓ<sup>2</sup>d) + ε<sub>−</sub><sup>6</sup>d<sup>3</sup>ℓ<sup>3</sup>k<sup>7</sup> m<sup>6</sup>ℓ<sup>12</sup> with constant

p

success probability of size

O(ε dℓk

−2

2

(m<sup>2</sup>ℓ<sup>2</sup>)<sup>1</sup><sub>/p</sub> log<sup>3</sup>(mℓ) log<sup>2</sup>(k) log(ε <sup>1</sup>) log(n)).

−

7

Conclusion

Our ﬁrst contribution involves investigating the VC dimension of range spaces characterized by

arbitrarily small perturbations of DTW distances. While our results hold for a relaxed variant

of the range spaces in question, they establish a robust link between numerous sampling

results dependent on the VC dimension and DTW distances. Indeed, our ﬁrst algorithmic

contribution is the construction of coresets for (k, ℓ)-median through the sensitivity sampling

framework by Feldman and Langberg [[25](#br25)]. Apart from the VC dimension, the crux of

adapting the sensitivity sampling framework to our (non-metric) setting was to use an

already known weak version of the triangle inequality satisﬁed by DTW. This inequality

prompted us to further explore approximation algorithms by approximating DTW with a

metric. By reducing to the metric case and plugging in our coresets, we designed an algorithm

for the (k, ℓ)-median problem, with running time linear in the number of the input sequences,

and an approximation factor predominantly determined by our generalised iterated triangle

inequality.

Although our primary motivation lies in constructing coresets, there are additional direct

consequences through sampling bounds that establish a connection between the sample size

and the VC dimension. For instance, suppose that we have a large set of time series, following

some unknown distribution, and we want to estimate the probability that a new time series

falls within a given DTW ball b. Suppose that we also allow for small perturbations of

the distances, i.e., we only want to guarantee that the estimated probability is realized by

some small perturbations of the distances. This probability can be approximated within

a constant additive error, by considering a random sample of size depending solely on the

VC dimension and the probability of success (over the random sampling) and measuring its

intersection with b (see e.g. Theorem [19).](#br11)[ ](#br11)Such an estimation can be used for example in

anomaly detection, where one aims to detect time series with a small chance of occurring, or

in time series segmentation, where diverse patterns may emerge throughout the series.

![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEACUDASIAAhEBAxEB/8QAGAABAAMBAAAAAAAAAAAAAAAAAAUGCAr/xAArEAAAAwYDCAMBAAAAAAAAAAABAgQAAwUGCBEHCRITFBUWGVmY1hchYXH/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A638QaBBxJnaZJ3GtXMAkMZuii2NnlPDao7laSZdO90nGGyrAeTV3B4QTVpdot7UaSlKG1G1xpEKyyQfQ5G9HMFzOyCdyU2h3VeJHZb3+iF5CGxfy4/1jGCQ6Yhe4RmfeWI+gs6Yhe4RmfeWI+gsYwaXptpzPT6ee4QOPdSeOJI8olxWRfUPij8lxOA7ijiLoyaXFfAoLwxIvFSD2Iudm/wB6fJkp9RNjYxjGD//Z)![ref9]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEACADASIAAhEBAxEB/8QAFgABAQEAAAAAAAAAAAAAAAAAAAYK/8QAKBAAAAMHBAICAwAAAAAAAAAAAQIGAAMEBQgREgcJFiEUFRMiUWFx/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/ANb+oNAg6krZSLca1dwBBirppGzs6T02qO4siU6d7icZalZDw2O9PKCZYu4Ly4jEpSh8o2uMRKtskH0ug3o7gu52QTuSmwd1XiR2W9+iF4ENi/q4/wBYxgTXbJB1Lop4G4LudnEjq+Dyq/MhvsULGKKCC4d/kO7dte6d0EDpgt06uwrSr81BMkpi4mhEhqfUZyxDKE5nT0nhqlP8OgPby8oGyLDeXD2MADn1ZjGD/9k=)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEACUDASIAAhEBAxEB/8QAGAABAAMBAAAAAAAAAAAAAAAAAAUGCAr/xAArEAAAAwYDCAMBAAAAAAAAAAABAgQAAwUGCBEHCRITFBUWGVmY1hchYXH/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A638QaBBxJnaZJ3GtXMAkMZuii2NnlPDao7laSZdO90nGGyrAeTV3B4QTVpdot7UaSlKG1G1xpEKyyQfQ5G9HMFzOyCdyU2h3VeJHZb3+iF5CGxfy4/1jGCQ6Yhe4RmfeWI+gs6Yhe4RmfeWI+gsYwaXptpzPT6ee4QOPdSeOJI8olxWRfUPij8lxOA7ijiLoyaXFfAoLwxIvFSD2Iudm/wB6fJkp9RNjYxjGD//Z)![ref5]![ref16]

<a name="br24"></a> 

24

Fast Approximations and Coresets for (k, ℓ)-Median under Dynamic Time Warping

References

1

Waleed H. Abdulla, David Chow, and Gary Sin. Cross-words reference template for dtw-based

speech recognition systems. In TENCON 2003. Conference on Convergent Technologies for

Asia-Paciﬁc Region, volume 4, pages 1576–1579 Vol.4, 2003.

2

3

4

Marcel R. Ackermann, Johannes Blömer, and Christian Sohler. Clustering for metric and

nonmetric distance measures. ACM Transactions on Algorithms, 6(4):59:1–59:26, 2010.

Martin Anthony and Peter L. Bartlett. Neural Network Learning: Theoretical Foundations.

Cambridge University Press, 1999. [doi:10.1017/CBO9780511624216](https://doi.org/10.1017/CBO9780511624216).

Vijay Arya, Naveen Garg, Rohit Khandekar, Adam Meyerson, Kamesh Munagala, and

Vinayaka Pandit. Local Search Heuristics for k-Median and Facility Location Problems. SIAM

Journal on Computing, 33(3):544–562, 2004.

5

6

Anselm Blumer, A. Ehrenfeucht, David Haussler, and Manfred K. Warmuth. Learnability and

the Vapnik-Chervonenkis dimension. Journal of the ACM, 36(4):929–965, oct 1989. URL:

<https://dl.acm.org/doi/10.1145/76359.76371>, [doi:10.1145/76359.76371](https://doi.org/10.1145/76359.76371).

Milutin Brankovic, Kevin Buchin, Koen Klaren, André Nusser, Aleksandr Popov, and

Sampson Wong. (k, l)-Medians Clustering of Trajectories Using Continuous Dynamic Time

Warping. In Proceedings of the 28th International Conference on Advances in Geographic

Information Systems, volume 1, pages 99–110, New York, NY, USA, nov 2020. ACM. URL:

<https://dl.acm.org/doi/10.1145/3397536.3422245>, [arXiv:arXiv:2012.00464v1](http://arxiv.org/abs/arXiv:2012.00464v1), [doi:10.](https://doi.org/10.1145/3397536.3422245)

[1145/3397536.3422245](https://doi.org/10.1145/3397536.3422245).

7

8

Vladimir Braverman, Vincent Cohen-Addad, Shaofeng H.-C. Jiang, Robert Krauthgamer,

Chris Schwiegelshohn, Mads Bech Toftrup, and Xuan Wu. The power of uniform sampling for

coresets. In 63rd IEEE Annual Symposium on Foundations of Computer Science, FOCS 2022,

Denver, CO, USA, October 31 - November 3, 2022, pages 462–473. IEEE, 2022.

Vladimir Braverman, Vincent Cohen-Addad, Shaofeng H.-C. Jiang, Robert Krauthgamer,

Chris Schwiegelshohn, Mads Bech Toftrup, and Xuan Wu. The power of uniform sampling

for coresets. In 63rd IEEE Annual Symposium on Foundations of Computer Science, FOCS

2022, Denver, CO, USA, October 31 - November 3, 2022, pages 462–473. IEEE, 2022. [doi:](https://doi.org/10.1109/FOCS54457.2022.00051)

[10.1109/FOCS54457.2022.00051](https://doi.org/10.1109/FOCS54457.2022.00051).

9

Vladimir Braverman, Dan Feldman, and Harry Lang. New Frameworks for Oﬄine and

Streaming Coreset Constructions. CoRR, abs/1612.00889, 2016. [arXiv:1612.00889](http://arxiv.org/abs/1612.00889).

10 Markus Brill, Till Fluschnik, Vincent Froese, Brijnesh J. Jain, Rolf Niedermeier, and David

Schultz. Exact mean computation in dynamic time warping spaces. Data Min. Knowl. Discov.,

33(1):252–291, 2019.

11 Kevin Buchin, Anne Driemel, Joachim Gudmundsson, Michael Horton, Irina Kostitsyna,

Maarten Löﬄer, and Martijn Struijs. Approximating (k, ℓ)-center clustering for curves. In

Timothy M. Chan, editor, Proceedings of the Thirtieth Annual ACM-SIAM Symposium on

Discrete Algorithms, SODA, pages 2922–2938, San Diego, California, USA, January 2019.

SIAM.

12 Kevin Buchin, Anne Driemel, and Martijn Struijs. On the hardness of computing an average

curve. In 17th Scandinavian Symposium and Workshops on Algorithm Theory, SWAT 2020,

June 22-24, 2020, Tórshavn, Faroe Islands, pages 19:1–19:19, 2020.

13 Kevin Buchin, Anne Driemel, and Martijn Struijs. On the Hardness of Computing an Average

Curve. In Susanne Albers, editor, 17th Scandinavian Symposium and Workshops on Algorithm

Theory, volume 162 of LIPIcs, pages 19:1–19:19, Tórshavn, Faroe Islands, June 2020. Schloss

Dagstuhl - Leibniz-Zentrum für Informatik.

14 Kevin Buchin, Anne Driemel, Natasja van de L’Isle, and André Nusser. klcluster: Center-

based Clustering of Trajectories. In Proceedings of the 27<sup>th</sup> ACM SIGSPATIAL International

Conference on Advances in Geographic Information Systems, pages 496–499, 2019.

15 Maike Buchin, Anne Driemel, and Dennis Rohde. Approximating (k, ℓ)-median clustering for

polygonal curves. In Dániel Marx, editor, Proceedings of the 2021 ACM-SIAM Symposium on

![](data:image/jpeg;base64,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)

<a name="br25"></a> 

Conradi, Kolbe, Psarros and Rohde

25

Discrete Algorithms, SODA 2021, Virtual Conference, January 10 - 13, 2021, pages 2697–2717.

SIAM, 2021.

16 Maike Buchin, Anne Driemel, and Dennis Rohde. Approximating (k, ℓ)-Median Clustering

for Polygonal Curves. In Dániel Marx, editor, Proceedings of the ACM-SIAM Symposium on

Discrete Algorithms, SODA, pages 2697–2717, Virtual Conference, January 2021. SIAM.

17 Maike Buchin, Anne Driemel, and Dennis Rohde. Approximating (k, ℓ)-median clustering for

polygonal curves. ACM Trans. Algorithms, 19(1):4:1–4:32, 2023.

18 Maike Buchin, Anne Driemel, Koen van Greevenbroek, Ioannis Psarros, and Dennis Rohde.

Approximating length-restricted means under dynamic time warping. In Parinya Chalermsook

and Bundit Laekhanukit, editors, Approximation and Online Algorithms - 20th International

Workshop, WAOA 2022, Potsdam, Germany, September 8-9, 2022, Proceedings, volume 13538

of Lecture Notes in Computer Science, pages 225–253. Springer, 2022.

19 Maike Buchin and Dennis Rohde. Coresets for (k, ℓ)-Median Clustering Under the Fréchet

Distance. In Niranjan Balachandran and R. Inkulu, editors, Algorithms and Discrete Applied

Mathematics - 8<sup>th</sup> International Conference, CALDAM, Puducherry, India, February 10-12,

Proceedings, volume 13179 of Lecture Notes in Computer Science, pages 167–180. Springer,

2022\.

20 Laurent Bulteau, Vincent Froese, and Rolf Niedermeier. Tight hardness results for consensus

problems on circular strings and time series. SIAM J. Discret. Math., 34(3):1854–1883, 2020.

21 Ke Chen. On Coresets for k-Median and k-Means Clustering in Metric and Euclidean Spaces

and Their Applications. SIAM Journal on Computing, 39(3):923–947, 2009.

22 Siu-Wing Cheng and Haoqiang Huang. Curve simpliﬁcation and clustering under fréchet

distance. In Nikhil Bansal and Viswanath Nagarajan, editors, Proceedings of the 2023 ACM-

SIAM Symposium on Discrete Algorithms, SODA 2023, Florence, Italy, January 22-25, 2023,

pages 1414–1432. SIAM, 2023.

23 Michael B. Cohen, Yin Tat Lee, Gary Miller, Jakub Pachocki, and Aaron Sidford. Geometric

median in nearly linear time. In Proceedings of the Forty-Eighth Annual ACM Symposium

on Theory of Computing, STOC ’16, page 9–21, New York, NY, USA, 2016. Association for

Computing Machinery. [doi:10.1145/2897518.2897647](https://doi.org/10.1145/2897518.2897647).

24 Anne Driemel, Amer Krivosija, and Christian Sohler. Clustering time series under the Fréchet

distance. In Proceedings of the Twenty-Seventh Annual ACM-SIAM Symposium on Discrete

Algorithms, pages 766–785, 2016.

25 Dan Feldman and Michael Langberg. A uniﬁed framework for approximating and clustering

data. In Lance Fortnow and Salil P. Vadhan, editors, Proceedings of the 43<sup>rd</sup> ACM Symposium

on Theory of Computing, pages 569–578. ACM, 2011.

26 Dan Feldman, Melanie Schmidt, and Christian Sohler. Turning Big Data Into Tiny Data:

Constant-Size Coresets for k-Means, PCA, and Projective Clustering. SIAM Journal on

Computing, 49(3):601–657, 2020.

27 Sariel Har-Peled and Micha Sharir. Relative (p,ε)-Approximations in Geometry. Discrete &

Computational Geometry, 45(3):462–496, April 2011.

28 Ville Hautamäki, Pekka Nykanen, and Pasi Franti. Time-series clustering by approximate

prototypes. In 2008 19th International Conference on Pattern Recognition, pages 1–4, 2008.

29 Piotr Indyk. Sublinear time algorithms for metric space problems. In Jeﬀrey Scott Vitter,

Lawrence L. Larmore, and Frank Thomson Leighton, editors, Proceedings of the Thirty-First

Annual ACM Symposium on Theory of Computing, May 1-4, 1999, Atlanta, Georgia, USA,

pages 428–434. ACM, 1999. [doi:10.1145/301250.301366](https://doi.org/10.1145/301250.301366).

30 Youngseon Jeong, Myong Kee Jeong, and Olufemi A. Omitaomu. Weighted dynamic time

warping for time series classiﬁcation. Pattern Recognit., 44(9):2231–2240, 2011.

31 Rohit J. Kate. Using dynamic time warping distances as features for improved time series

classiﬁcation. Data Min. Knowl. Discov., 30(2):283–312, 2016.

32 Amit Kumar, Yogish Sabharwal, and Sandeep Sen. A Simple Linear Time (1+ε)-Approximation

Algorithm for k-Means Clustering in Any Dimensions. In 45th Symposium on Foundations of

![ref16]

<a name="br26"></a> 

26

Fast Approximations and Coresets for (k, ℓ)-Median under Dynamic Time Warping

Computer Science (FOCS), 17-19 October, Rome, Italy, Proceedings, pages 454–462. IEEE

Computer Society, 2004.

33 Michael Langberg and Leonard J. Schulman. Universal ε-approximators for Integrals. In

Proceedings of the 21<sup>st</sup> Annual ACM-SIAM Symposium on Discrete Algorithms (SODA), pages

598–607, 2010.

34 Daniel Lemire. Faster retrieval with a two-pass dynamic-time-warping lower bound. Pattern

Recognition, 42(9):2169 – 2180, 2009.

35 Alexander Munteanu, Chris Schwiegelshohn, Christian Sohler, and David P. Woodruﬀ. On

Coresets for Logistic Regression. In Samy Bengio, Hanna M. Wallach, Hugo Larochelle,

Kristen Grauman, Nicolò Cesa-Bianchi, and Roman Garnett, editors, Advances in Neural

Information Processing Systems 31: Annual Conference on Neural Information Processing

Systems, NeurIPS, December 3-8, Montréal, Canada, pages 6562–6571, 2018.

36 François Petitjean, Germain Forestier, Geoﬀrey I. Webb, Ann E. Nicholson, Yanping Chen,

and Eamonn J. Keogh. Dynamic time warping averaging of time series allows faster and more

accurate classiﬁcation. In Ravi Kumar, Hannu Toivonen, Jian Pei, Joshua Zhexue Huang,

and Xindong Wu, editors, 2014 IEEE International Conference on Data Mining, ICDM 2014,

Shenzhen, China, December 14-17, 2014, pages 470–479. IEEE Computer Society, 2014.

37 François Petitjean, Germain Forestier, Geoﬀrey I. Webb, Ann E. Nicholson, Yanping Chen,

and Eamonn J. Keogh. Faster and more accurate classiﬁcation of time series by exploiting a

novel dynamic time warping averaging algorithm. Knowl. Inf. Syst., 47(1):1–26, 2016.

38 François Petitjean, Alain Ketterlin, and Pierre Gançarski. A global averaging method for

dynamic time warping, with applications to clustering. Pattern Recognit., 44(3):678–693, 2011.

39 Lawrence Rabiner and Jay Wilpon. Considerations in applying clustering techniques to speaker

independent word recognition. In ICASSP ’79. IEEE International Conference on Acoustics,

Speech, and Signal Processing, volume 4, pages 578–581, 1979.

40 Thanawin Rakthanmanon, Bilson J. L. Campana, Abdullah Mueen, Gustavo E. A. P. A.

Batista, M. Brandon Westover, Qiang Zhu, Jesin Zakaria, and Eamonn J. Keogh. Addressing

big data time series: Mining trillions of time series subsequences under dynamic time warping.

ACM Trans. Knowl. Discov. Data, 7(3):10:1–10:31, 2013.

41 Norbert Sauer. On the density of families of sets. Journal of Combinatorial Theory Series A,

13:145–147, 1972.

42 Nathan Schaar, Vincent Froese, and Rolf Niedermeier. Faster binary mean computation under

dynamic time warping. In 31st Annual Symposium on Combinatorial Pattern Matching, CPM

2020, June 17-19, 2020, Copenhagen, Denmark, pages 28:1–28:13, 2020.

43 Saharon Shelah. A combinatorial problem; stability and order for models and theories in

inﬁnitary languages. Paciﬁc Journal of Mathematics, 41(1), 1972.

44 Tuan Minh Tran, Xuan-May Thi Le, Hien T. Nguyen, and Van-Nam Huynh. A novel non-

parametric method for time series classiﬁcation based on k-nearest neighbors and dynamic

time warping barycenter averaging. Eng. Appl. Artif. Intell., 78:173–185, 2019.

45 Vladimir Vapnik and Alexey Chervonenkis. On the uniform convergence of relative frequencies

of events to their probabilities. Theory of Probability and its Applications, 16:264–280, 1971.

A

VC dimension analysis

ꢀ

ꢋ

ꢌꢁ

▶ Lemma 43. If p is even, then the VC dimension of X<sup>d</sup> , B<sub>r</sub><sup>p</sup><sub>,m</sub>(σ) | σ ∈ X<sup>d</sup> is in

=m

\=

ℓ

ꢀ

ꢁ

O dℓ mp

<sup>2</sup> log( ) .

Proof. For any two vectors x = (x , . . . , x ), y = (y , . . . , y ) let x ⊕ y ∈ R<sup>2</sup> be their

d

1

d

1

d

concatenation (x , . . . , x , y , . . . , y ). We apply Theorem [7](#br6)[ ](#br6)as follows. The parameter vector

1

d

1

d

α(σ, r) encodes the center σ and the radius r of each ball. More formally, for any ball

center σ = (σ , . . . , σ ) and radius r, we deﬁne a polynomial function parameterized by

1

ℓ

α(σ, r) = σ<sub>1</sub> ⊕ · · · ⊕ σ ⊕ (r) ∈ R<sup>dℓ</sup>

<sup>+1</sup>. Let be the family of at most m<sup>O(ℓ) polynomial</sup>

F

ℓ

![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAWABwDASIAAhEBAxEB/8QAGwAAAQQDAAAAAAAAAAAAAAAAAAEGBwkCCAr/xAAxEAABAgQFAwICCwAAAAAAAAABAhEDBAUGAAcSITEIIkFRYRMjFBUWMjNCUoGRobH/xAAXAQADAQAAAAAAAAAAAAAAAAAGBwgK/8QAJhEAAQMDAwMFAQAAAAAAAAAAAQIDBAUGEQAHEhMhMQgUIkFhFv/aAAwDAQACEQMRAD8A6iczszr/AKbmHdUjJ3fXabS6dWFplZNCwpKklS3Qn5qWQNIbtffgcYZULNrMkICze9wQzEPxdOrfvHP47hz7HbGObaQjMy8vzfCq8+UlbrYwlpCfvOS2tT+u3oMR3pZSgD58g7ewc8M+wO38vli3Q3O3Cjbg3eP7C6wlN2TAEioAAAYx2Ch9fgxn9xq/LUta2H7ZgvPU+M665GjrccXHbUvmpptSu6kE+cq7qPk4ye+pFi5s5kkdt9V/W/a6wQo/pV8/YE+XPHGLBen646xVstKZPVycnZ+oxZ6qCNMxoo1rCJtaUA9xYABgHLeNmxVnytmDFwdvQeCzvx6NxizLpkgoj5S0hcRKVK+sqwlyNyEzqwOCPGKc9Cd/7h1/dC4WZFbmTYrNmOrbeqcnrSOoKnS0qTy5EcQCcY8DQNvNblEpds0+ZBiMMdSpx2T0mm21HMaSohRQEEjKM+SAe+o9vjpXr9zXbXrgl7ro8tAq85NTUKXjyc6uNBTMqSpKIi0HQpSQliU7HxhrHo7uYkn7ZUNj6yFQ/d+4AP7f0d8GDFM3j6cNlKlcdZmTrFiPyZVWkSZDpqtwILj6yOThS1VkIST9hKUp/NLulbgXfDpkeLGrCm2G2Gkob9lTlBKUNpCRyXDUs4AHlRJx3zpE9HVyhQJvGhkAnb6BPs2lhs5/3wMbcZQ5dTeXtkSNsT1RlqjMSs3PzCpqUgxIUFaZuYMZKQiL3hSArSoklyHG2DBhgbDbJbX2XdFQnWxazVJlPURyK66zVa49zYMmI4UFEqpvtj5oQrkEBfxxyxkEevS8bkrlOjRarU1y47cht1DRjQ2kpcShaAoFiO0o4StQwSR3zjODr//Z)

<a name="br27"></a> 

Conradi, Kolbe, Psarros and Rohde

27

Figure 7 Illustration of a lopsided traversal.

functions, each one realizing the diﬀerence between the p<sup>th</sup> power of the radius and the p<sup>th</sup>

power of the cost of a diﬀerent traversal. More formally, for a curve τ = (τ , . . . , τ ) ∈ X :

d

m

1

m









X

F = α(σ, r), τ<sub>1</sub> ⊕ · · · ⊕ τ → r<sup>p</sup> −

∥σ − τ ∥ | σ ∈ , r ∈

p

X<sup>d R</sup>+, T ∈ T

.

m

i

j

2

ℓ

<sup>ℓ,m</sup>



(i,j)∈T

Since p is even, the square root in ∥ · ∥ vanishes and the degree is upper bounded by p.

2

Finally, we deﬁne a boolean function g : {−1, 1} → {0, 1} which simply returns 1 iﬀ

|F |

at least one of the arguments is non-negative. The function g deﬁnes H = (X, R<sup>p</sup> ) as a

m,ℓ

<sup>(</sup><sub>ℓ</sub><sup>)</sup>-combination of sign( ), and by Theorem [7,](#br6)[ ](#br6)we conclude that the VC dimension of

m<sup>O</sup>

F

H

◀

ꢀ

ꢁ

are in O dℓ<sup>2</sup> log(mp) .

B

Approximating the DTW distance by a metric

▶ Corollary 44. Let X be a set of curves of complexity at most m. Let k and ℓ be given. There

/p

is an algorithm which computes a ((4mℓ)<sup>1</sup> (41 +O(ε)), 1)-approximation to the (k, ℓ)-median

problem on X under dtw in O(nm<sup>3</sup>ℓ log<sup>3</sup>(mε<sub>−</sub><sup>1</sup>) + n<sup>2</sup>ℓ<sup>2</sup> + n<sup>3</sup> + nk + k<sup>7</sup>ε<sub>−</sub><sup>5</sup> log<sup>5</sup> n).

Proof. This is a direct consequence of Theorem [29,](#br18)[ ](#br18)Lemma [26](#br17)[ ](#br17)and Theorem [49.](#br28)

◀

▶ Corollary 45. Let X be a set of curves of complexity at most m. Let k and ℓ be given. There

/p

is an algorithm which computes a ((4mℓ)<sup>1</sup> (62 +O(ε)), 1)-approximation to the (k, ℓ)-median

problem on X under dtw<sub>p</sub> in O(nm<sup>3</sup> + n<sup>2</sup>ℓ<sup>2</sup> + n<sup>3</sup> + nk + k<sup>7</sup>ε<sub>−</sub><sup>5</sup> log<sup>5</sup> n).

Proof. This is a direct consequence of Theorem [29,](#br18)[ ](#br18)Lemma [26](#br17)[ ](#br17)and Proposition [37.](#br21)

◀

C

Approximate ℓ-simpliﬁcations under DTW

▶ Deﬁnition 46. With the notation of Deﬁnition [1,](#br4)[ ](#br4)we cal l a traversal lopsided (refer to

Figure [7)](#br27)[ ](#br27)if for all i we have that (a<sub>i+1</sub>, b<sub>i+1</sub>) − (a , b ) ∈ {(0, 1), (1, 1)}.

i

i

▶ Lemma 47. Let σ ∈ X<sup>d</sup> and τ ∈ X<sup>d</sup> . There is an optimal traversal T∗ ∈ T that realizes

m

n

n,m

dtw(σ, τ) and for all 1 ≤ i ≤ |T | − 2 it holds that (a , b ) − (a , b )̸= (1, 1).

∗

i+2 i+2

i

i

Proof. Let T be an arbitrary traversal for which there is an index i such that (a<sub>i+2</sub>, b ) −

i+2

(a , b ) = (1, 1). Observe that we can remove the index pair (a<sub>i+1</sub>, b<sub>i+1</sub>) from the traversal

i

i

without increasing the distance, which implies the claim.

◀

![](data:image/jpeg;base64,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)![ref3]![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAWABwDASIAAhEBAxEB/8QAGwAAAQQDAAAAAAAAAAAAAAAAAAEGBwkCCAr/xAAxEAABAgQFAwICCwAAAAAAAAABAhEDBAUGAAcSITEIIkFRYRMjFBUWMjNCUoGRobH/xAAXAQADAQAAAAAAAAAAAAAAAAAGBwgK/8QAJhEAAQMDAwMFAQAAAAAAAAAAAQIDBAUGEQAHEhMhMQgUIkFhFv/aAAwDAQACEQMRAD8A6iczszr/AKbmHdUjJ3fXabS6dWFplZNCwpKklS3Qn5qWQNIbtffgcYZULNrMkICze9wQzEPxdOrfvHP47hz7HbGObaQjMy8vzfCq8+UlbrYwlpCfvOS2tT+u3oMR3pZSgD58g7ewc8M+wO38vli3Q3O3Cjbg3eP7C6wlN2TAEioAAAYx2Ch9fgxn9xq/LUta2H7ZgvPU+M665GjrccXHbUvmpptSu6kE+cq7qPk4ye+pFi5s5kkdt9V/W/a6wQo/pV8/YE+XPHGLBen646xVstKZPVycnZ+oxZ6qCNMxoo1rCJtaUA9xYABgHLeNmxVnytmDFwdvQeCzvx6NxizLpkgoj5S0hcRKVK+sqwlyNyEzqwOCPGKc9Cd/7h1/dC4WZFbmTYrNmOrbeqcnrSOoKnS0qTy5EcQCcY8DQNvNblEpds0+ZBiMMdSpx2T0mm21HMaSohRQEEjKM+SAe+o9vjpXr9zXbXrgl7ro8tAq85NTUKXjyc6uNBTMqSpKIi0HQpSQliU7HxhrHo7uYkn7ZUNj6yFQ/d+4AP7f0d8GDFM3j6cNlKlcdZmTrFiPyZVWkSZDpqtwILj6yOThS1VkIST9hKUp/NLulbgXfDpkeLGrCm2G2Gkob9lTlBKUNpCRyXDUs4AHlRJx3zpE9HVyhQJvGhkAnb6BPs2lhs5/3wMbcZQ5dTeXtkSNsT1RlqjMSs3PzCpqUgxIUFaZuYMZKQiL3hSArSoklyHG2DBhgbDbJbX2XdFQnWxazVJlPURyK66zVa49zYMmI4UFEqpvtj5oQrkEBfxxyxkEevS8bkrlOjRarU1y47cht1DRjQ2kpcShaAoFiO0o4StQwSR3zjODr//Z)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAWABsDASIAAhEBAxEB/8QAGwAAAQQDAAAAAAAAAAAAAAAAAAEGBwkCCAr/xAAvEAABAgQFAQcDBQAAAAAAAAABAgQDBQYRAAcSITFBCBMiI0JRcRYkgTJTYZHB/8QAFwEAAwEAAAAAAAAAAAAAAAAABgcICv/EACoRAAICAgEDAgUFAQAAAAAAAAECAwQFBhEAByESEwgUMUFRIjJhcYGh/9oADAMBAAIRAxEAPwDqJzNzOr+W5h1UyZ1fPpdLJdOYiGjKGvUlQUIvlpPeosgaRbbrxfYMqFmzmSEBX1vUEPvPM0hdrd4b/ve5J/j43xjm2kIzMrL193N3xQV+OyocRCUfqve3eKBJ52PQWju2kqF9gevT8/ji+wAxlh7n9zu4VbuDt4G4bYEXbb4AXIgAceoccerj+hx4/wCdX7quq6zPrOPmmx9aSSSrXkkkkrxM5dokZvLIfuSfLH6/f6GRYubOZRA0V1P9fpSpVwr3B8/w7dSDY9NsWQ5NT6ZTHLKk3s1cu3kwcMnS3TmJFBXFWJk9SFKOo7hCUp54AxUudyRYE6VHpa4QSCDa/Pz/AJi1XIaBDjZR0VEWlJUpg8ube02mA6EdAMVl8BO+dw8/ve7R2s1buVoNQqPBPkbHu2Gc5agrgtyQQoPAI/aBx/PS872a9hcVg8Pao1YYBPkhEwihiiY8VLDcN6ApIBT7kgf4eobrfsrz+pqtn9QQKrlDaBN3jp1CgRmb1caClwuGtKVrQdClICDcp2JtbjDVPY7qUlR+spHY8fYPyed+ttwNxa1/jcwYoPcPhw7K5LYc1du6LTns2ctaszynK7Ahkndj6pCseWRFJ/Cqq/gdAmM7g7hTxlarWzDR144YlSP5LGsFVUQKAz02bwFH1Yk8eeT0iex1UoJJrGRnwrTuxf8AqQUpPPIJuf75xublrRDiiqHp+l3b6A9cShtHgxXTaHEhwIxjPXToGGiL40hKY6UHVuVJJGxGDBhmdgey/bPSc9mrOr6xHiJreGjrWHhymbmEkK3K8qoUtZOdF4kRW9SKrEjy3BII3u+3bFnadOvlsk1uGCz7sUZr04gkntvH6h8vXiJ/QzDgkjzzxz56/9k=)![ref2]

<a name="br28"></a> 

28

Fast Approximations and Coresets for (k, ℓ)-Median under Dynamic Time Warping

Figure 8 Illustration of Proof of Lemma [48.](#br27)[ ](#br27)The vertex that can safely be deleted is marked on

the left, and removed on the right.

▶ Lemma 48. Let τ ∈ X<sup>d</sup> and σ ∈ X<sup>d</sup> satisfy

m

ℓ

dtw (σ, τ) = inf dtw (σ , τ).

(1)

′

p

p

σ<sup>′</sup>∈X<sup>d</sup>

ℓ

Then dtw<sub>p</sub>(σ, τ) can be realized with a lopsided traversal.

Proof. Refer to Figure [8.](#br28)[ ](#br28)Assume that dtw (σ, τ) is not attained by a lopsided traversal. It

p

suﬃces to show that there is then a curve σ ∈ X with dtw (σ , τ) ≤ dtw(σ, τ).

′

d

′

p

ℓ

Fix some traversal. There is an index i such that (a , b ) − (a , b ) = (1, 0). Observe

i+1 i+1

i

i

that removing vertex σ<sub>a</sub> from σ does not increase dtw (σ, τ). Indeed, by Lemma [47,](#br27)[ ](#br27)we

p

i

+1

can assume that (a<sub>i+2</sub>, b ) − (a , b )̸= (0, 1), so a = a + 1. Therefore, discarding

i+2

i+1 i+1

i+2

i+1

σ<sub>a</sub> removes exactly one summand from the sum in dtw(σ, τ), yielding a contradiction.

◀

i

+1

By the following theorem, we can calculate a (1 + ε)-approximate ℓ-simpliﬁcation for any ε

for the classical DTW distance dtw := dtw<sub>1</sub> with some success probability.

▶ Theorem 49. There is a constant δ < 1 σ ∈ X<sup>d</sup> , ℓ ∈ N and ε > 0, one can compute in

m

O(ℓm<sup>2</sup> + m<sup>3</sup>d log<sup>3</sup>(mε−<sup>1</sup>)) time a curve

<sub>d</sub> such that

ℓ

∗

σ ∈ X

inf dtw(σ , σ) ≤ dtw(σ , σ) ≤ (1 + ε) inf dtw(σ , σ),

∗

ℓ

ℓ

σ ∈X<sup>d</sup>

σ ∈X<sup>d</sup>

ℓ

ℓ

ℓ

ℓ

with success probability (1 − δ)<sub>ℓ</sub>.

Proof. If m ≤ ℓ, return σ. If m > ℓ, the idea is to use Lemma [48](#br27)[ ](#br27)and ﬁnd a curve in X

d

ℓ

along with a lopsided traversal such that the inequality holds.

If σ ∈ X is such that dtw(σ , σ) is minimal over all curves in X , by Lemma [48,](#br27)[ ](#br27)there

∗

d

∗

d

ℓ

ℓ

exists a lopsided traversal realizing the distance. Observe that the lopsided traversal induces

a partitioning of the vertices of σ into ℓ distinct sets G , . . . , G of vertices. Each G contains

1

ℓ

i

a list of consecutive vertices of σ that connect to a given vertex of σ = (σ , . . . , σ ) and we

∗

∗

1

∗

ℓ

have that

X

ℓ

X

dtw(σ<sub>∗</sub>, σ) =

∥σ<sub>i</sub> − p∥.

i=1 p∈G

i

Conversely, any partitioning of the vertices of σ into ℓ distinct sets of consecutive vertices

together with a point in R for each of these sets gives rise to a curve in X . We thus see

d

d

ℓ

that it is enough to ﬁnd such a partition along with a vertex for every set in the partition,

such that the inequality holds.

This problem can be solved similarly to one dimensional ℓ-medians clustering. That is, for

P

σ = (σ1, . . . , σ ), we deﬁne a cost function C (a, b) = inf

d

b

∥σ − c∥

, corresponding

m

σ

c∈R i=a

i

to the cost of grouping indices σ , . . . , σ , with some optimal point in R . It is well-

d

a

b

known that this cost function can be (1 + ε)-approximated in O(md log3(mε<sub>−</sub>1)) time [[23](#br25)]

![](data:image/jpeg;base64,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)![](data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAANADIDASIAAhEBAxEB/8QAGwAAAQQDAAAAAAAAAAAAAAAACQACBgsECAr/xAAoEAABBAICAQMEAwEAAAAAAAABAgMEBQYHCBESAAkTFBYhMRgiQYH/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A7yLW2q6hl6Zf2kasr2mZSn1WUppmChCkhY+odeU222oJacUyx/ZXSXfEno+hqP8AuYY/sdUmm4a6d2xzDtZiBTYxtXBsfcp+Jy8/jqcRKxzOd4pVdLwBNOkoVkk1OG5Aa4S64CPI+oPhEKTBcX9xzde/qvkpTw8t0Nxz3blmm8H48zW1y8Gu8z169EXYbjz51TjByjIJabCvRg1QK6ub10G8i8bfKPuYipLbU0dLXV7ESrqoVZCS0UNw4MdqNHbSUoQoJaZQhHakNNpUop7UEJB/Q9AK1HIH3fvFCVe3ZxoK0qAcV/PKbFaW4O/6sj+NbvytH8kK7T8v4PSOvWbE5fc99cvIyLkr7dsyLrdavoo6+JO5pPKDaqsmkBTlYiZq9zWWr0RsQcisWSrvJvuhxVVOTURRWSxZl2KVtEZltAbbQEISPFKUgAJB676H+d9Dv0xcVkoKQCgBtSAEnoAKIJ/HX77SP+E+gHZr33NuKWbZ1jOnMxy234773yiemFRcd+QlMNc7nsWn470irtY2KOTbNpmmvIkeTLorFdkoT4kd9z4WSjxO8LeTUT0mLDi3FJKfKhHaiotYTkh1opUsmMG3nFzCktp+RpLbZ/PyeY8OjEdoab1fvHEcv1rtzBMYz7CcmrYFXkmP5DUx5cK7rXFmQa+d2EvrhpkRo7yGUPp8VMpAV++6733yuAuR+2Ds2u9w/izvCu1NYUm9anFtWab1Hqp7WmL67rocLJXoUo3TWysjXkFtKg1Yg5LPVQ04yBcl6X9PXkhghZIiwjEAmYkE/khLpIBP5IB6HYH6B6HY/wAHpeg/cUua26dy8XONm38vexxOW7W0Fp3ZOUJqahcCqTkWc68x3KLtNZBVMfVDrxZ2koQopfeMeN8bJdcKPMr0H//Z)![](data:image/jpeg;base64,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)![ref3]

<a name="br29"></a> 

Conradi, Kolbe, Psarros and Rohde

29

with constant success probability 1 − δ. Hence, it can be approximately computed in

O(m d mε

<sup>3</sup> log<sup>3</sup>( −<sup>1</sup>)) time for all index pairs and , such that the algorithm succeeds in

a

b

approximating the cost for an optimal partition with probability (1 − δ) , where δ < 1 is the

ℓ

constant success probability of the algorithm presented in [[23](#br25)]. We next compute D(i, j), the

approximate minimal cost to group the vertices σ , . . . , σ into j groups for all i ≤ m, j ≤ ℓ.

1

i

This can be done in O(ℓm<sup>2</sup>) time, by ﬁrst setting D(i, 1) = C(1, i) and then computing

D(i, j) = min D(l, j − 1) + C(l + 1, i). Overall, it thus takes O(ℓm<sup>2 +</sup> m d mε<sup>−1))</sup>

<sup>3</sup> log<sup>3</sup>(

l≤i

time to (1 + ε)-approximately compute an optimal curve with success probability 1 − δ .

◀

ℓ

▶ Corollary 50. For σ ∈ X<sup>d</sup> and m ≥ ℓ > 0, one can compute in O(ℓm<sup>2</sup> + m<sup>3</sup>) time a curve

m

with vertices a subset of vertices in σ, such that

′

σ

dtw(σ , σ) = inf dtw(σ , σ),

′

ℓ

σ

ℓ

where the inﬁmum is taken over all σ ∈ X with vertices lying on vertices of σ.

d

ℓ

ℓ

Proof. This follows immediately from the fact that the cost C (a, b) of clustering the a-th to

P

b-th vertex of σ with a single vertex of σ can be computed in O(m(b−a)) = O(m<sup>2) time.</sup>

◀

C.1 Approximate ℓ-simpliﬁcations under dtw<sub>p</sub>

The same ideas as above lead to a deterministic 2-approximation for dtw for arbitrary values

p

of p.

▶ Proposition 37. For σ = (σ , . . . , σ ) ∈ X<sup>d</sup> and integer ℓ > 0, one can compute in

1

m

<sub>d</sub> such that

ℓ

m

O(m d

<sup>2</sup>( + + )) time a curve

ℓ

m

σ ∈ X

∗

inf dtw (σ , σ) ≤ dtw (σ , σ) ≤ 2 inf dtw (σ , σ).

∗

p

ℓ

p

p

ℓ

σ ∈X<sup>d</sup>

σ ∈X<sup>d</sup>

ℓ

ℓ

ℓ

ℓ

Proof. First compute in O(m<sup>2</sup>d) time all pairwise distances of vertices of σ, and store them

P

for later use. Next compute C(a, b) = min<sub>a≤i≤b</sub>

∥σ − σ ∥<sup>p</sup> for all 1 ≤ a ≤ b ≤ m.

a≤j≤b

i

j

2

This can be done in O(m<sup>3</sup>) time, as for any ﬁxed 1 ≤ a ≤ b ≤ m, the computations of C(a, b)

are a big part of the computations of C(a, b + 1). To be more precise, when computing

P

C(a, b) we store O(m) temporary values ν(i) =

C(a, b) = min

∥σ − σ ∥

<sub>p</sub> for all

a ≤ i ≤ b

. Then

i

j

2

a≤j≤b

ν(i). Once we have these values, C(a, b + 1) can be computed in O(m)

a≤i≤b

time, by updating ν(i) with ν(i) + ∥p − p ∥<sup>p</sup>. We then compute ν(b + 1) in O(m) time.

i

b+1

2

As C(a, a) can be computed in O(1) time for any 1 ≤ a ≤ m, we can compute all values

C(a, b) for 1 ≤ a ≤ b ≤ m in O(m

<sup>3</sup>) time. Consider, for 1

≤ ℓ ≤ m ≤ m

,

′

′

′

X<sup>ℓ</sup>

D(ℓ<sup>′</sup>, m<sub>′</sub>) =

min

1

C(s<sub>i−1</sub> + 1, s<sub>i</sub>),

<sub>σ</sub>′=(σ ,...,σ

′

′

)

′

=1

m

i

S=(s1,...,s ) subsequence of σ

′

′

ℓ

\=

σ

m

s

′

′

ℓ

with s := 0, which corresponds to the optimal partitioning of the ﬁrst m indices into

′

0

ℓ

contiguous disjoint intervals under the cost function C. Computing the value D(ℓ, m)

′

takes O(m<sup>2</sup>ℓ) time via the recursive formula D(ℓ , m ) = min D(ℓ − 1, l) + C(l + 1, m ).

′

′

′

′

l≤m′

The subsequence realizing D(ℓ, m) deﬁnes σ . For correctness, observe that the optimal

∗

simpliﬁcation also yields a partition of the vertices of σ into ℓ contiguous disjoint intervals.

Let T be the partition computed, and let T<sup>opt</sup> = ([a , b ], . . . , [a , b ]) be an optimal partition

1

1

ℓ

ℓ

![ref2]

<a name="br30"></a> 

30

Fast Approximations and Coresets for (k, ℓ)-Median under Dynamic Time Warping

of [m], that realizes inf

d

dtw (σ , σ). For a vertex v in a ﬁxed optimal simpliﬁcation,

p

ℓ

i

σ ∈X

ℓ

let π be a closest point amo<sup>ℓ</sup>ng {σ , . . . , σ }. Then

i

a

b

i

i









1/p

1/p

X

X

X



C a, b

(

) = 

min

a≤i≤b

<sup>p</sup>

∥σ − σ ∥

i

j

2

[a,b]∈T

[a,b]∈T

X

a≤j≤b





1/p

X

≤ 

min

a≤i≤b

∥σ − σ ∥

<sup>p</sup>

i

j

2

[a,b]∈T <sup>opt</sup>

X

a≤j≤b





1/p

X

a≤j≤b

X

≤ 

∥π − σ ∥

<sup>p</sup>

i

j

2

[a,b]∈T <sup>opt</sup>

X





1/p

≤ 

(∥π − v ∥<sub>2</sub> + ∥v − σ ∥<sub>2</sub>)

<sup>p</sup>

i

i

i

j

[a,b]∈T <sup>opt</sup>

X

a≤j≤b

X





1/p

≤ 

(2∥v − σ ∥<sub>2</sub>)

<sup>p</sup>

i

j

[a,b]∈T <sup>opt</sup>

a≤j≤b





1/p

X

X

≤ 2 

∥v − σ ∥

<sup>p</sup> = 2 inf dtw (

)

◀

i

j

2

p σ , σ .

ℓ

σ ∈X<sup>d</sup>

[a,b]∈T <sup>opt</sup>

ℓ

a≤j≤b

ℓ



[ref1]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAPAA8DASIAAhEBAxEB/8QAFgABAQEAAAAAAAAAAAAAAAAABQYK/8QAIxAAAQQCAQUAAwAAAAAAAAAAAQIDBAUGBxIIERMUIQAVIv/EABUBAQEAAAAAAAAAAAAAAAAAAAYI/8QAIBEAAgICAgIDAAAAAAAAAAAAAQIDBAYRACESFTFBUf/aAAwDAQACEQMRAD8A2M7O23sl3ay9Ua0exuFaMYO9maLa9aXMgWL6FpCYy4SFxzAb+EqnexJI5fI/8nvLdP299lZ7nsjGM2r8dnxZzF+/jmVYvzNQt3HlQUXdc28pJ9zw/sK8omAsl4LIMdHEEs7a1FtdG3W9qatmY0LSTjwxeY1kQUptuu5KKkxEoUPE+sLPN76FcUHgOP4L046I2dg2YKuMws6OLEp6ibX1FRTLVKjl+4W0uzuZzhDPntLIwoglvBtHL12+4Pb5M9i1lL5wtWxNlscU2VXIbVKpYj9J6VI42oeBKggkGQkBR3ro64caSzaso7K4VWBDaAH1rokgDxDH5Hej+8//2Q==
[ref2]: data:image/jpeg;base64,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
[ref3]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAMAAwDASIAAhEBAxEB/8QAFgABAQEAAAAAAAAAAAAAAAAABwQK/8QAIRAAAgIDAAICAwAAAAAAAAAAAQIDBAUGEQchABITFDH/xAAVAQEBAAAAAAAAAAAAAAAAAAAICf/EACQRAAEDAgQHAAAAAAAAAAAAAAIBBAUDBgcSIjEAERMUMmGh/9oADAMBAAIRAxEAPwDUH5O8nb/R8g7XVrbbsVDHY/OWo6GPitFkKn8wVVCzjirwcHAB6/h+P3gncNmzGn3rOZzd27aj2C3AkttmkmSA43E2EiZw7AhZLErLw+g3Dwg/LNr8Aabm9mzWctZPZord+d7Ukde5iv145pGIdoUsYWxKv2DMOPNJwH1w+/iDoXjnB6lhpsXi7OTeu9+S2725qc0zTS1KcbkulGFefWFAFCADh564BEiyrXxMjMVbzn5+9Ul4yTR8DCOB3IITIkeNyQiztqYeIqOglXb3wip6Vt53bESwj4btHLXoLXdLTboVYkpcj1ipmSEa5tfxd//Z
[ref4]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAHAAwDASIAAhEBAxEB/8QAFgABAQEAAAAAAAAAAAAAAAAAAAkK/8QAIBAAAQQDAAIDAAAAAAAAAAAABAECAwUABgcIESEiJP/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDSz2vnvn5ddc6BZc52fcxtGM2wmPU68LplHTV4tVEsqxfgm2EciBzkexFY4ZrkRPlqesoD4i0nZKDjgQHeCrQzoS396QbPb7ELtBSgTzQuAa21DNsIXQNYkiQjoQroG/VzGe0TGMD/2Q==
[ref5]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAWABkDASIAAhEBAxEB/8QAGwAAAQQDAAAAAAAAAAAAAAAAAAYHCAkCBQr/xAAzEAABAgUDAQYBDQAAAAAAAAABAhEDBAUGIQAHEjEIExZBUWHRFCMyM0JSU3GRoaKx4f/EABkBAAIDAQAAAAAAAAAAAAAAAAYIBQcJCv/EACURAAEEAQMEAgMAAAAAAAAAAAIBAwQFBgcRIQASMUEIIhMUUf/aAAwDAQACEQMRAD8A6itzdzr/AKbuFdUlJ3fXqbTKdWIiZWTQsKSpKjFBQg96kBA4hsDqceWkTC3Z3JCAvxvcEMxXi8QtvrC+fnnyXf0/LOsd2wIe5l5fa7urz3ErdbKhxEBH0n6c1OSQ+D5Bm8KeKlhwwJ92z0f0YdPIAa5YtT9TtQo2oOXomYZWIjlk5BFLEdkROOPttwi+NuOeOUTp/MVxXGX8Zr3nq+M465Gjmbhx2iNTJoCLciBfaqvKr59+3coW5+4M7cFvy0xe1XjS0eu0eFMwZlPeJjQIlRlkqhJaOSgqcAqbCeQbJ1bPzh/iK/l8NUrWsCbotsl3FwUUukkAkVKVGfVurH39H1dE6/f9P81oL8B87vbXEM7k2k+daSUvqgP2bF788hQSsJRDvRV+o+BT1x/V6pXXKrra2wx1mtjtRmSr5ZkLLbbXeavtIpH2CiEW2ybrvsnG/HUML47K9fua7a9cEvddHloFXnJqahy8eSnVxoKZhaFJStaCEKUnhnjglmZtJY9ju5ipR8ZUPJJzIT5OfXIGQA+G/vRo0X5l8cNFLLIrqbNwWG/Jk2smS+6trkAK4+a/ZxRbtgBFXfwIiP8AB6jKvUHL4dZGixrgm47bLQg3+lWmgiIAgohHDI+ERPJKq7c9bKidkW46dWaTUIl30SKiQqchOqQmRnwuIiUmoMwpCSVcQpYhlIJDAlyWfU+vkJ++P3+GjRphPjvpRp9hNJfw8WxxipjTLGG/Jabm2khHXW434wNSmTpJiogqpsBCK+VRV56B82yS7vpUB22nFLOPGJplVYjM9gE42RCgx2WRXdfZIq+kXbjr/9k=
[ref6]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEAAQDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAUEAEAAAAAAAAAAAAAAAAAAAAA/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AL+AAf/Z
[ref7]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEAFEDASIAAhEBAxEB/8QAGAABAQADAAAAAAAAAAAAAAAAAAUGCQr/xAAuEAAABAIKAQIHAQAAAAAAAAAAAQIDBAkFBggREhMZWZjWFgcUFSEyQVFmleT/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A63/UKwIfqTXSs1dztrTAahqrbSUdTa6p+m1o7xapNXlumlR0bVagfDY74RRDePC3Be7iMKUpLNO688IoqWST1GwTpzBZnaDch21YG7V5obReX0oT4EeFJfYrzu/IAAoaYidwiZ9yxPoIaYidwiZ9yxPoIAAaYidwiZ9yxPoIaYidwiZ9yxPoIAAaYidwiZ9yxPoIaYidwiZ9yxPoIAAaYidwiZ9yxPoIhRktz2kapgre0yiI9siFiidirUua8+p1x0ktRLngqc1iHyTOHbuLLU8+eI8z5AAbFfCP2+vH9/8AyAAAP//Z
[ref8]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEABYDASIAAhEBAxEB/8QAFwABAAMAAAAAAAAAAAAAAAAAAAYICv/EACcQAAAEBAUEAwAAAAAAAAAAAAECBQYDBAcSAAgJERQTFRZBFxgi/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/ANb9Qsgg1Jejme451tQFhmdqlPLZ2nTbMd4syW8eKJTCmtZB8Nnu0JEO+2HJcuYtKUodUdtxhKVpkBFTpKKOoNqdkGJLQjWQ819kMu5AG0hQYX5KHoNx2D3hhgLD0AydFoE51hzfaTONWrvCCdC8fr/W/wCR2wl3qEgod5R0kWwj8JeJweERR5EW1PnJ+W6I8i8jDDAf/9k=
[ref9]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEAAkDASIAAhEBAxEB/8QAFgABAQEAAAAAAAAAAAAAAAAAAAgK/8QAJRAAAAQEBQUAAAAAAAAAAAAAAgMFBgAEBwgBERITFhQVVZTj/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/ANcVULC+dVCcT9LvPv5ZAnCsKDmEz6fXFcbYqMMZhJpaIgN7h0525vSuQi5VM6w7ZLMMDvi1Z4Uton/Ornv/AChCA//Z
[ref10]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEAAsDASIAAhEBAxEB/8QAFgABAQEAAAAAAAAAAAAAAAAAAAgK/8QAIxAAAQMCBgMBAAAAAAAAAAAABgECBAUHAAgSExQVAxEWIf/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDWtcbIU+5BsZH7s6uf8GcSzvJX2h1ucx3zAGOPnPc51JFh742d1FEi6ESHA5kjYark3n+/yxLfDkm2YSMADS0yPWidIi0VDK5VcUpPSNIjNPalZFxoHcVmT71S53DjbzkRdpuGGA//2Q==
[ref11]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEAC8DASIAAhEBAxEB/8QAGAABAAMBAAAAAAAAAAAAAAAAAAQGBwr/xAApEAAAAwYFBAMBAAAAAAAAAAABAgYAAwQFBxIICREUFhlZmNYTFUEi/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/AOt+oWAQakrRTLcca2YCgzK2ZR07Ok6bYjuLIlPHeiUwy1LSHhsd9RKHd9ruC3cRaUpQ+UdNRpMqyyAey6CejmDZnZBeQzo1jvFfY7LqQBtIUEF/JQ/A1HQP1jGCf0xC9wjM+8sR9BZ0xC9wjM+8sR9BYxgdMQvcIzPvLEfQW0ijuCc1D6opZblxd44KtkhYdRQh0ZW+vvPkFG7+TvoQsTHJ3ikp3EXLzPd5LX28JtY126f2nssExg//2Q==
[ref12]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEADUDASIAAhEBAxEB/8QAGAABAQADAAAAAAAAAAAAAAAAAAUGCAr/xAAtEAAAAwYDCAEFAAAAAAAAAAABAgMABAUGBwkIERITFRYZQVmY1hQXISMyUf/EABQBAQAAAAAAAAAAAAAAAAAAAAD/xAAUEQEAAAAAAAAAAAAAAAAAAAAA/9oADAMBAAIRAxEAPwDrfqFgEGpM6TNO441rgMhmm2JP0bPKdNsR3C0ky8dUSmGGytAeDX7dEIT16U3L5bxpKUobUcsxwiFWyQWhrkqNwW52QVHdM2hPFeJEyZh+pC8BDpKHQMxy/rGMFDliF7hFz7yxH0FnLEL3CLn3liPoLGMDliF7hFz7yxH0Fpx7ZYEiQlC4Jc6MBIeZcBPiu1CYSrgGzMPAX3SNlmYnUeoMYwbWUIoOrh5l+MS6jW7EFWjf0WJHTx2v9SPqPMcLMLk7uW6oLFNzQb4MFArsDyDhsFcntZdba/k0lMYwf//Z
[ref13]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEAAgDASIAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAr/xAAiEAAABAYCAwAAAAAAAAAAAAACAwUGAAEEBwgSExUUFhf/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8ArwemEIla7KjdwrMDOZDGN2SfvzZAv/11pgGEq9KoBaQGb6kdqxxgCNMMQOxnsknH0fly35JIQgP/2Q==
[ref14]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEAC0DASIAAhEBAxEB/8QAGAABAAMBAAAAAAAAAAAAAAAAAAUGCAr/xAArEAAAAwUGBQUAAAAAAAAAAAABAgMABQYHEQQICRIVFhQZWZjWExchYXH/xAAUAQEAAAAAAAAAAAAAAAAAAAAA/8QAFBEBAAAAAAAAAAAAAAAAAAAAAP/aAAwDAQACEQMRAD8A635g3BBmTG0SRuN9XEAgMYueltfZ4TlteO2tBMOnVynF2wq4dm27R3QTNlTsXF2jKUpQ9UaVGkOrDJBZ3WNUcQXE7IJ0SmyJ3rxImWtfghdhDQv1Uf1jGCQ5YheoRifd2I+As5YheoRifd2I+AsYwOWIXqEYn3diPgLaEkbdGLI9GJUgvM3uZvDEqrpUFaeU5/cBdw6SR4kBOG1NtunS0njqAmepKL8YaxO8ap8NQ5jB/9k=
[ref15]: data:image/jpeg;base64,/9j/4AAQSkZJRgABAQEAYABgAAD/2wBDAAEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/2wBDAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQH/wAARCAAEABcDASIAAhEBAxEB/8QAFwABAAMAAAAAAAAAAAAAAAAAAAYHCv/EACcQAAAEBAUFAQEAAAAAAAAAAAECBQYDBBESAAcICRQTFRYiQRgk/8QAFAEBAAAAAAAAAAAAAAAAAAAAAP/EABQRAQAAAAAAAAAAAAAAAAAAAAD/2gAMAwEAAhEDEQA/ANb2YWgQcyXm53wOtXcAYYu1Rnlw7Sy21HeLMhvHi2mFNayB4bPdoSCX2w5LlzFpSlDqjSowpK2yAip0lFHcG3OyDEloRrIeq+yGWpAG0hQYXqUPgVGgfcMMBY2XO30XLp5ozy/bO4M+xRu40a2Y+pLyhmqYqCVPJVVdDFmSXNGSCeFQT/6oXGVJSSm/fj9MzDDAf//Z
[ref16]: data:image/jpeg;base64,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
