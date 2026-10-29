/*
 * Copyright (c) 2026 harmor.
 *
 * AES-128/192/256 encryption core and ECB/CBC wrappers.
 *
 * The block function dispatches to the RISC-V Zvkned implementation when
 * AES_USE_ZVKNED is enabled.  The portable C path is constant time: it uses
 * a bitsliced S-box built from the GF(2^8) inverse x^254 and has no secret-
 * dependent loads or branches.
 */

#include "aes_ct.h"

#include <string.h>

static const uint8_t aes_rcon[11] = {
    0x00, 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36
};

static const uint8_t aes_sbox_tt[256] = {
    0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5,
    0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76,
    0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0,
    0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0,
    0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc,
    0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15,
    0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a,
    0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75,
    0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0,
    0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84,
    0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b,
    0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf,
    0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85,
    0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8,
    0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5,
    0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2,
    0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17,
    0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73,
    0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88,
    0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb,
    0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c,
    0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79,
    0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9,
    0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08,
    0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6,
    0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a,
    0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e,
    0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e,
    0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94,
    0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf,
    0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68,
    0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16,
};

static const uint8_t aes_inv_sbox_tt[256] = {
    0x52, 0x09, 0x6a, 0xd5, 0x30, 0x36, 0xa5, 0x38,
    0xbf, 0x40, 0xa3, 0x9e, 0x81, 0xf3, 0xd7, 0xfb,
    0x7c, 0xe3, 0x39, 0x82, 0x9b, 0x2f, 0xff, 0x87,
    0x34, 0x8e, 0x43, 0x44, 0xc4, 0xde, 0xe9, 0xcb,
    0x54, 0x7b, 0x94, 0x32, 0xa6, 0xc2, 0x23, 0x3d,
    0xee, 0x4c, 0x95, 0x0b, 0x42, 0xfa, 0xc3, 0x4e,
    0x08, 0x2e, 0xa1, 0x66, 0x28, 0xd9, 0x24, 0xb2,
    0x76, 0x5b, 0xa2, 0x49, 0x6d, 0x8b, 0xd1, 0x25,
    0x72, 0xf8, 0xf6, 0x64, 0x86, 0x68, 0x98, 0x16,
    0xd4, 0xa4, 0x5c, 0xcc, 0x5d, 0x65, 0xb6, 0x92,
    0x6c, 0x70, 0x48, 0x50, 0xfd, 0xed, 0xb9, 0xda,
    0x5e, 0x15, 0x46, 0x57, 0xa7, 0x8d, 0x9d, 0x84,
    0x90, 0xd8, 0xab, 0x00, 0x8c, 0xbc, 0xd3, 0x0a,
    0xf7, 0xe4, 0x58, 0x05, 0xb8, 0xb3, 0x45, 0x06,
    0xd0, 0x2c, 0x1e, 0x8f, 0xca, 0x3f, 0x0f, 0x02,
    0xc1, 0xaf, 0xbd, 0x03, 0x01, 0x13, 0x8a, 0x6b,
    0x3a, 0x91, 0x11, 0x41, 0x4f, 0x67, 0xdc, 0xea,
    0x97, 0xf2, 0xcf, 0xce, 0xf0, 0xb4, 0xe6, 0x73,
    0x96, 0xac, 0x74, 0x22, 0xe7, 0xad, 0x35, 0x85,
    0xe2, 0xf9, 0x37, 0xe8, 0x1c, 0x75, 0xdf, 0x6e,
    0x47, 0xf1, 0x1a, 0x71, 0x1d, 0x29, 0xc5, 0x89,
    0x6f, 0xb7, 0x62, 0x0e, 0xaa, 0x18, 0xbe, 0x1b,
    0xfc, 0x56, 0x3e, 0x4b, 0xc6, 0xd2, 0x79, 0x20,
    0x9a, 0xdb, 0xc0, 0xfe, 0x78, 0xcd, 0x5a, 0xf4,
    0x1f, 0xdd, 0xa8, 0x33, 0x88, 0x07, 0xc7, 0x31,
    0xb1, 0x12, 0x10, 0x59, 0x27, 0x80, 0xec, 0x5f,
    0x60, 0x51, 0x7f, 0xa9, 0x19, 0xb5, 0x4a, 0x0d,
    0x2d, 0xe5, 0x7a, 0x9f, 0x93, 0xc9, 0x9c, 0xef,
    0xa0, 0xe0, 0x3b, 0x4d, 0xae, 0x2a, 0xf5, 0xb0,
    0xc8, 0xeb, 0xbb, 0x3c, 0x83, 0x53, 0x99, 0x61,
    0x17, 0x2b, 0x04, 0x7e, 0xba, 0x77, 0xd6, 0x26,
    0xe1, 0x69, 0x14, 0x63, 0x55, 0x21, 0x0c, 0x7d,
};

static const uint32_t aes_te0[256] = {
    0xc66363a5u, 0xf87c7c84u, 0xee777799u, 0xf67b7b8du,
    0xfff2f20du, 0xd66b6bbdu, 0xde6f6fb1u, 0x91c5c554u,
    0x60303050u, 0x02010103u, 0xce6767a9u, 0x562b2b7du,
    0xe7fefe19u, 0xb5d7d762u, 0x4dababe6u, 0xec76769au,
    0x8fcaca45u, 0x1f82829du, 0x89c9c940u, 0xfa7d7d87u,
    0xeffafa15u, 0xb25959ebu, 0x8e4747c9u, 0xfbf0f00bu,
    0x41adadecu, 0xb3d4d467u, 0x5fa2a2fdu, 0x45afafeau,
    0x239c9cbfu, 0x53a4a4f7u, 0xe4727296u, 0x9bc0c05bu,
    0x75b7b7c2u, 0xe1fdfd1cu, 0x3d9393aeu, 0x4c26266au,
    0x6c36365au, 0x7e3f3f41u, 0xf5f7f702u, 0x83cccc4fu,
    0x6834345cu, 0x51a5a5f4u, 0xd1e5e534u, 0xf9f1f108u,
    0xe2717193u, 0xabd8d873u, 0x62313153u, 0x2a15153fu,
    0x0804040cu, 0x95c7c752u, 0x46232365u, 0x9dc3c35eu,
    0x30181828u, 0x379696a1u, 0x0a05050fu, 0x2f9a9ab5u,
    0x0e070709u, 0x24121236u, 0x1b80809bu, 0xdfe2e23du,
    0xcdebeb26u, 0x4e272769u, 0x7fb2b2cdu, 0xea75759fu,
    0x1209091bu, 0x1d83839eu, 0x582c2c74u, 0x341a1a2eu,
    0x361b1b2du, 0xdc6e6eb2u, 0xb45a5aeeu, 0x5ba0a0fbu,
    0xa45252f6u, 0x763b3b4du, 0xb7d6d661u, 0x7db3b3ceu,
    0x5229297bu, 0xdde3e33eu, 0x5e2f2f71u, 0x13848497u,
    0xa65353f5u, 0xb9d1d168u, 0x00000000u, 0xc1eded2cu,
    0x40202060u, 0xe3fcfc1fu, 0x79b1b1c8u, 0xb65b5bedu,
    0xd46a6abeu, 0x8dcbcb46u, 0x67bebed9u, 0x7239394bu,
    0x944a4adeu, 0x984c4cd4u, 0xb05858e8u, 0x85cfcf4au,
    0xbbd0d06bu, 0xc5efef2au, 0x4faaaae5u, 0xedfbfb16u,
    0x864343c5u, 0x9a4d4dd7u, 0x66333355u, 0x11858594u,
    0x8a4545cfu, 0xe9f9f910u, 0x04020206u, 0xfe7f7f81u,
    0xa05050f0u, 0x783c3c44u, 0x259f9fbau, 0x4ba8a8e3u,
    0xa25151f3u, 0x5da3a3feu, 0x804040c0u, 0x058f8f8au,
    0x3f9292adu, 0x219d9dbcu, 0x70383848u, 0xf1f5f504u,
    0x63bcbcdfu, 0x77b6b6c1u, 0xafdada75u, 0x42212163u,
    0x20101030u, 0xe5ffff1au, 0xfdf3f30eu, 0xbfd2d26du,
    0x81cdcd4cu, 0x180c0c14u, 0x26131335u, 0xc3ecec2fu,
    0xbe5f5fe1u, 0x359797a2u, 0x884444ccu, 0x2e171739u,
    0x93c4c457u, 0x55a7a7f2u, 0xfc7e7e82u, 0x7a3d3d47u,
    0xc86464acu, 0xba5d5de7u, 0x3219192bu, 0xe6737395u,
    0xc06060a0u, 0x19818198u, 0x9e4f4fd1u, 0xa3dcdc7fu,
    0x44222266u, 0x542a2a7eu, 0x3b9090abu, 0x0b888883u,
    0x8c4646cau, 0xc7eeee29u, 0x6bb8b8d3u, 0x2814143cu,
    0xa7dede79u, 0xbc5e5ee2u, 0x160b0b1du, 0xaddbdb76u,
    0xdbe0e03bu, 0x64323256u, 0x743a3a4eu, 0x140a0a1eu,
    0x924949dbu, 0x0c06060au, 0x4824246cu, 0xb85c5ce4u,
    0x9fc2c25du, 0xbdd3d36eu, 0x43acacefu, 0xc46262a6u,
    0x399191a8u, 0x319595a4u, 0xd3e4e437u, 0xf279798bu,
    0xd5e7e732u, 0x8bc8c843u, 0x6e373759u, 0xda6d6db7u,
    0x018d8d8cu, 0xb1d5d564u, 0x9c4e4ed2u, 0x49a9a9e0u,
    0xd86c6cb4u, 0xac5656fau, 0xf3f4f407u, 0xcfeaea25u,
    0xca6565afu, 0xf47a7a8eu, 0x47aeaee9u, 0x10080818u,
    0x6fbabad5u, 0xf0787888u, 0x4a25256fu, 0x5c2e2e72u,
    0x381c1c24u, 0x57a6a6f1u, 0x73b4b4c7u, 0x97c6c651u,
    0xcbe8e823u, 0xa1dddd7cu, 0xe874749cu, 0x3e1f1f21u,
    0x964b4bddu, 0x61bdbddcu, 0x0d8b8b86u, 0x0f8a8a85u,
    0xe0707090u, 0x7c3e3e42u, 0x71b5b5c4u, 0xcc6666aau,
    0x904848d8u, 0x06030305u, 0xf7f6f601u, 0x1c0e0e12u,
    0xc26161a3u, 0x6a35355fu, 0xae5757f9u, 0x69b9b9d0u,
    0x17868691u, 0x99c1c158u, 0x3a1d1d27u, 0x279e9eb9u,
    0xd9e1e138u, 0xebf8f813u, 0x2b9898b3u, 0x22111133u,
    0xd26969bbu, 0xa9d9d970u, 0x078e8e89u, 0x339494a7u,
    0x2d9b9bb6u, 0x3c1e1e22u, 0x15878792u, 0xc9e9e920u,
    0x87cece49u, 0xaa5555ffu, 0x50282878u, 0xa5dfdf7au,
    0x038c8c8fu, 0x59a1a1f8u, 0x09898980u, 0x1a0d0d17u,
    0x65bfbfdau, 0xd7e6e631u, 0x844242c6u, 0xd06868b8u,
    0x824141c3u, 0x299999b0u, 0x5a2d2d77u, 0x1e0f0f11u,
    0x7bb0b0cbu, 0xa85454fcu, 0x6dbbbbd6u, 0x2c16163au,
};

static const uint32_t aes_te1[256] = {
    0xa5c66363u, 0x84f87c7cu, 0x99ee7777u, 0x8df67b7bu,
    0x0dfff2f2u, 0xbdd66b6bu, 0xb1de6f6fu, 0x5491c5c5u,
    0x50603030u, 0x03020101u, 0xa9ce6767u, 0x7d562b2bu,
    0x19e7fefeu, 0x62b5d7d7u, 0xe64dababu, 0x9aec7676u,
    0x458fcacau, 0x9d1f8282u, 0x4089c9c9u, 0x87fa7d7du,
    0x15effafau, 0xebb25959u, 0xc98e4747u, 0x0bfbf0f0u,
    0xec41adadu, 0x67b3d4d4u, 0xfd5fa2a2u, 0xea45afafu,
    0xbf239c9cu, 0xf753a4a4u, 0x96e47272u, 0x5b9bc0c0u,
    0xc275b7b7u, 0x1ce1fdfdu, 0xae3d9393u, 0x6a4c2626u,
    0x5a6c3636u, 0x417e3f3fu, 0x02f5f7f7u, 0x4f83ccccu,
    0x5c683434u, 0xf451a5a5u, 0x34d1e5e5u, 0x08f9f1f1u,
    0x93e27171u, 0x73abd8d8u, 0x53623131u, 0x3f2a1515u,
    0x0c080404u, 0x5295c7c7u, 0x65462323u, 0x5e9dc3c3u,
    0x28301818u, 0xa1379696u, 0x0f0a0505u, 0xb52f9a9au,
    0x090e0707u, 0x36241212u, 0x9b1b8080u, 0x3ddfe2e2u,
    0x26cdebebu, 0x694e2727u, 0xcd7fb2b2u, 0x9fea7575u,
    0x1b120909u, 0x9e1d8383u, 0x74582c2cu, 0x2e341a1au,
    0x2d361b1bu, 0xb2dc6e6eu, 0xeeb45a5au, 0xfb5ba0a0u,
    0xf6a45252u, 0x4d763b3bu, 0x61b7d6d6u, 0xce7db3b3u,
    0x7b522929u, 0x3edde3e3u, 0x715e2f2fu, 0x97138484u,
    0xf5a65353u, 0x68b9d1d1u, 0x00000000u, 0x2cc1ededu,
    0x60402020u, 0x1fe3fcfcu, 0xc879b1b1u, 0xedb65b5bu,
    0xbed46a6au, 0x468dcbcbu, 0xd967bebeu, 0x4b723939u,
    0xde944a4au, 0xd4984c4cu, 0xe8b05858u, 0x4a85cfcfu,
    0x6bbbd0d0u, 0x2ac5efefu, 0xe54faaaau, 0x16edfbfbu,
    0xc5864343u, 0xd79a4d4du, 0x55663333u, 0x94118585u,
    0xcf8a4545u, 0x10e9f9f9u, 0x06040202u, 0x81fe7f7fu,
    0xf0a05050u, 0x44783c3cu, 0xba259f9fu, 0xe34ba8a8u,
    0xf3a25151u, 0xfe5da3a3u, 0xc0804040u, 0x8a058f8fu,
    0xad3f9292u, 0xbc219d9du, 0x48703838u, 0x04f1f5f5u,
    0xdf63bcbcu, 0xc177b6b6u, 0x75afdadau, 0x63422121u,
    0x30201010u, 0x1ae5ffffu, 0x0efdf3f3u, 0x6dbfd2d2u,
    0x4c81cdcdu, 0x14180c0cu, 0x35261313u, 0x2fc3ececu,
    0xe1be5f5fu, 0xa2359797u, 0xcc884444u, 0x392e1717u,
    0x5793c4c4u, 0xf255a7a7u, 0x82fc7e7eu, 0x477a3d3du,
    0xacc86464u, 0xe7ba5d5du, 0x2b321919u, 0x95e67373u,
    0xa0c06060u, 0x98198181u, 0xd19e4f4fu, 0x7fa3dcdcu,
    0x66442222u, 0x7e542a2au, 0xab3b9090u, 0x830b8888u,
    0xca8c4646u, 0x29c7eeeeu, 0xd36bb8b8u, 0x3c281414u,
    0x79a7dedeu, 0xe2bc5e5eu, 0x1d160b0bu, 0x76addbdbu,
    0x3bdbe0e0u, 0x56643232u, 0x4e743a3au, 0x1e140a0au,
    0xdb924949u, 0x0a0c0606u, 0x6c482424u, 0xe4b85c5cu,
    0x5d9fc2c2u, 0x6ebdd3d3u, 0xef43acacu, 0xa6c46262u,
    0xa8399191u, 0xa4319595u, 0x37d3e4e4u, 0x8bf27979u,
    0x32d5e7e7u, 0x438bc8c8u, 0x596e3737u, 0xb7da6d6du,
    0x8c018d8du, 0x64b1d5d5u, 0xd29c4e4eu, 0xe049a9a9u,
    0xb4d86c6cu, 0xfaac5656u, 0x07f3f4f4u, 0x25cfeaeau,
    0xafca6565u, 0x8ef47a7au, 0xe947aeaeu, 0x18100808u,
    0xd56fbabau, 0x88f07878u, 0x6f4a2525u, 0x725c2e2eu,
    0x24381c1cu, 0xf157a6a6u, 0xc773b4b4u, 0x5197c6c6u,
    0x23cbe8e8u, 0x7ca1ddddu, 0x9ce87474u, 0x213e1f1fu,
    0xdd964b4bu, 0xdc61bdbdu, 0x860d8b8bu, 0x850f8a8au,
    0x90e07070u, 0x427c3e3eu, 0xc471b5b5u, 0xaacc6666u,
    0xd8904848u, 0x05060303u, 0x01f7f6f6u, 0x121c0e0eu,
    0xa3c26161u, 0x5f6a3535u, 0xf9ae5757u, 0xd069b9b9u,
    0x91178686u, 0x5899c1c1u, 0x273a1d1du, 0xb9279e9eu,
    0x38d9e1e1u, 0x13ebf8f8u, 0xb32b9898u, 0x33221111u,
    0xbbd26969u, 0x70a9d9d9u, 0x89078e8eu, 0xa7339494u,
    0xb62d9b9bu, 0x223c1e1eu, 0x92158787u, 0x20c9e9e9u,
    0x4987ceceu, 0xffaa5555u, 0x78502828u, 0x7aa5dfdfu,
    0x8f038c8cu, 0xf859a1a1u, 0x80098989u, 0x171a0d0du,
    0xda65bfbfu, 0x31d7e6e6u, 0xc6844242u, 0xb8d06868u,
    0xc3824141u, 0xb0299999u, 0x775a2d2du, 0x111e0f0fu,
    0xcb7bb0b0u, 0xfca85454u, 0xd66dbbbbu, 0x3a2c1616u,
};

static const uint32_t aes_te2[256] = {
    0x63a5c663u, 0x7c84f87cu, 0x7799ee77u, 0x7b8df67bu,
    0xf20dfff2u, 0x6bbdd66bu, 0x6fb1de6fu, 0xc55491c5u,
    0x30506030u, 0x01030201u, 0x67a9ce67u, 0x2b7d562bu,
    0xfe19e7feu, 0xd762b5d7u, 0xabe64dabu, 0x769aec76u,
    0xca458fcau, 0x829d1f82u, 0xc94089c9u, 0x7d87fa7du,
    0xfa15effau, 0x59ebb259u, 0x47c98e47u, 0xf00bfbf0u,
    0xadec41adu, 0xd467b3d4u, 0xa2fd5fa2u, 0xafea45afu,
    0x9cbf239cu, 0xa4f753a4u, 0x7296e472u, 0xc05b9bc0u,
    0xb7c275b7u, 0xfd1ce1fdu, 0x93ae3d93u, 0x266a4c26u,
    0x365a6c36u, 0x3f417e3fu, 0xf702f5f7u, 0xcc4f83ccu,
    0x345c6834u, 0xa5f451a5u, 0xe534d1e5u, 0xf108f9f1u,
    0x7193e271u, 0xd873abd8u, 0x31536231u, 0x153f2a15u,
    0x040c0804u, 0xc75295c7u, 0x23654623u, 0xc35e9dc3u,
    0x18283018u, 0x96a13796u, 0x050f0a05u, 0x9ab52f9au,
    0x07090e07u, 0x12362412u, 0x809b1b80u, 0xe23ddfe2u,
    0xeb26cdebu, 0x27694e27u, 0xb2cd7fb2u, 0x759fea75u,
    0x091b1209u, 0x839e1d83u, 0x2c74582cu, 0x1a2e341au,
    0x1b2d361bu, 0x6eb2dc6eu, 0x5aeeb45au, 0xa0fb5ba0u,
    0x52f6a452u, 0x3b4d763bu, 0xd661b7d6u, 0xb3ce7db3u,
    0x297b5229u, 0xe33edde3u, 0x2f715e2fu, 0x84971384u,
    0x53f5a653u, 0xd168b9d1u, 0x00000000u, 0xed2cc1edu,
    0x20604020u, 0xfc1fe3fcu, 0xb1c879b1u, 0x5bedb65bu,
    0x6abed46au, 0xcb468dcbu, 0xbed967beu, 0x394b7239u,
    0x4ade944au, 0x4cd4984cu, 0x58e8b058u, 0xcf4a85cfu,
    0xd06bbbd0u, 0xef2ac5efu, 0xaae54faau, 0xfb16edfbu,
    0x43c58643u, 0x4dd79a4du, 0x33556633u, 0x85941185u,
    0x45cf8a45u, 0xf910e9f9u, 0x02060402u, 0x7f81fe7fu,
    0x50f0a050u, 0x3c44783cu, 0x9fba259fu, 0xa8e34ba8u,
    0x51f3a251u, 0xa3fe5da3u, 0x40c08040u, 0x8f8a058fu,
    0x92ad3f92u, 0x9dbc219du, 0x38487038u, 0xf504f1f5u,
    0xbcdf63bcu, 0xb6c177b6u, 0xda75afdau, 0x21634221u,
    0x10302010u, 0xff1ae5ffu, 0xf30efdf3u, 0xd26dbfd2u,
    0xcd4c81cdu, 0x0c14180cu, 0x13352613u, 0xec2fc3ecu,
    0x5fe1be5fu, 0x97a23597u, 0x44cc8844u, 0x17392e17u,
    0xc45793c4u, 0xa7f255a7u, 0x7e82fc7eu, 0x3d477a3du,
    0x64acc864u, 0x5de7ba5du, 0x192b3219u, 0x7395e673u,
    0x60a0c060u, 0x81981981u, 0x4fd19e4fu, 0xdc7fa3dcu,
    0x22664422u, 0x2a7e542au, 0x90ab3b90u, 0x88830b88u,
    0x46ca8c46u, 0xee29c7eeu, 0xb8d36bb8u, 0x143c2814u,
    0xde79a7deu, 0x5ee2bc5eu, 0x0b1d160bu, 0xdb76addbu,
    0xe03bdbe0u, 0x32566432u, 0x3a4e743au, 0x0a1e140au,
    0x49db9249u, 0x060a0c06u, 0x246c4824u, 0x5ce4b85cu,
    0xc25d9fc2u, 0xd36ebdd3u, 0xacef43acu, 0x62a6c462u,
    0x91a83991u, 0x95a43195u, 0xe437d3e4u, 0x798bf279u,
    0xe732d5e7u, 0xc8438bc8u, 0x37596e37u, 0x6db7da6du,
    0x8d8c018du, 0xd564b1d5u, 0x4ed29c4eu, 0xa9e049a9u,
    0x6cb4d86cu, 0x56faac56u, 0xf407f3f4u, 0xea25cfeau,
    0x65afca65u, 0x7a8ef47au, 0xaee947aeu, 0x08181008u,
    0xbad56fbau, 0x7888f078u, 0x256f4a25u, 0x2e725c2eu,
    0x1c24381cu, 0xa6f157a6u, 0xb4c773b4u, 0xc65197c6u,
    0xe823cbe8u, 0xdd7ca1ddu, 0x749ce874u, 0x1f213e1fu,
    0x4bdd964bu, 0xbddc61bdu, 0x8b860d8bu, 0x8a850f8au,
    0x7090e070u, 0x3e427c3eu, 0xb5c471b5u, 0x66aacc66u,
    0x48d89048u, 0x03050603u, 0xf601f7f6u, 0x0e121c0eu,
    0x61a3c261u, 0x355f6a35u, 0x57f9ae57u, 0xb9d069b9u,
    0x86911786u, 0xc15899c1u, 0x1d273a1du, 0x9eb9279eu,
    0xe138d9e1u, 0xf813ebf8u, 0x98b32b98u, 0x11332211u,
    0x69bbd269u, 0xd970a9d9u, 0x8e89078eu, 0x94a73394u,
    0x9bb62d9bu, 0x1e223c1eu, 0x87921587u, 0xe920c9e9u,
    0xce4987ceu, 0x55ffaa55u, 0x28785028u, 0xdf7aa5dfu,
    0x8c8f038cu, 0xa1f859a1u, 0x89800989u, 0x0d171a0du,
    0xbfda65bfu, 0xe631d7e6u, 0x42c68442u, 0x68b8d068u,
    0x41c38241u, 0x99b02999u, 0x2d775a2du, 0x0f111e0fu,
    0xb0cb7bb0u, 0x54fca854u, 0xbbd66dbbu, 0x163a2c16u,
};

static const uint32_t aes_te3[256] = {
    0x6363a5c6u, 0x7c7c84f8u, 0x777799eeu, 0x7b7b8df6u,
    0xf2f20dffu, 0x6b6bbdd6u, 0x6f6fb1deu, 0xc5c55491u,
    0x30305060u, 0x01010302u, 0x6767a9ceu, 0x2b2b7d56u,
    0xfefe19e7u, 0xd7d762b5u, 0xababe64du, 0x76769aecu,
    0xcaca458fu, 0x82829d1fu, 0xc9c94089u, 0x7d7d87fau,
    0xfafa15efu, 0x5959ebb2u, 0x4747c98eu, 0xf0f00bfbu,
    0xadadec41u, 0xd4d467b3u, 0xa2a2fd5fu, 0xafafea45u,
    0x9c9cbf23u, 0xa4a4f753u, 0x727296e4u, 0xc0c05b9bu,
    0xb7b7c275u, 0xfdfd1ce1u, 0x9393ae3du, 0x26266a4cu,
    0x36365a6cu, 0x3f3f417eu, 0xf7f702f5u, 0xcccc4f83u,
    0x34345c68u, 0xa5a5f451u, 0xe5e534d1u, 0xf1f108f9u,
    0x717193e2u, 0xd8d873abu, 0x31315362u, 0x15153f2au,
    0x04040c08u, 0xc7c75295u, 0x23236546u, 0xc3c35e9du,
    0x18182830u, 0x9696a137u, 0x05050f0au, 0x9a9ab52fu,
    0x0707090eu, 0x12123624u, 0x80809b1bu, 0xe2e23ddfu,
    0xebeb26cdu, 0x2727694eu, 0xb2b2cd7fu, 0x75759feau,
    0x09091b12u, 0x83839e1du, 0x2c2c7458u, 0x1a1a2e34u,
    0x1b1b2d36u, 0x6e6eb2dcu, 0x5a5aeeb4u, 0xa0a0fb5bu,
    0x5252f6a4u, 0x3b3b4d76u, 0xd6d661b7u, 0xb3b3ce7du,
    0x29297b52u, 0xe3e33eddu, 0x2f2f715eu, 0x84849713u,
    0x5353f5a6u, 0xd1d168b9u, 0x00000000u, 0xeded2cc1u,
    0x20206040u, 0xfcfc1fe3u, 0xb1b1c879u, 0x5b5bedb6u,
    0x6a6abed4u, 0xcbcb468du, 0xbebed967u, 0x39394b72u,
    0x4a4ade94u, 0x4c4cd498u, 0x5858e8b0u, 0xcfcf4a85u,
    0xd0d06bbbu, 0xefef2ac5u, 0xaaaae54fu, 0xfbfb16edu,
    0x4343c586u, 0x4d4dd79au, 0x33335566u, 0x85859411u,
    0x4545cf8au, 0xf9f910e9u, 0x02020604u, 0x7f7f81feu,
    0x5050f0a0u, 0x3c3c4478u, 0x9f9fba25u, 0xa8a8e34bu,
    0x5151f3a2u, 0xa3a3fe5du, 0x4040c080u, 0x8f8f8a05u,
    0x9292ad3fu, 0x9d9dbc21u, 0x38384870u, 0xf5f504f1u,
    0xbcbcdf63u, 0xb6b6c177u, 0xdada75afu, 0x21216342u,
    0x10103020u, 0xffff1ae5u, 0xf3f30efdu, 0xd2d26dbfu,
    0xcdcd4c81u, 0x0c0c1418u, 0x13133526u, 0xecec2fc3u,
    0x5f5fe1beu, 0x9797a235u, 0x4444cc88u, 0x1717392eu,
    0xc4c45793u, 0xa7a7f255u, 0x7e7e82fcu, 0x3d3d477au,
    0x6464acc8u, 0x5d5de7bau, 0x19192b32u, 0x737395e6u,
    0x6060a0c0u, 0x81819819u, 0x4f4fd19eu, 0xdcdc7fa3u,
    0x22226644u, 0x2a2a7e54u, 0x9090ab3bu, 0x8888830bu,
    0x4646ca8cu, 0xeeee29c7u, 0xb8b8d36bu, 0x14143c28u,
    0xdede79a7u, 0x5e5ee2bcu, 0x0b0b1d16u, 0xdbdb76adu,
    0xe0e03bdbu, 0x32325664u, 0x3a3a4e74u, 0x0a0a1e14u,
    0x4949db92u, 0x06060a0cu, 0x24246c48u, 0x5c5ce4b8u,
    0xc2c25d9fu, 0xd3d36ebdu, 0xacacef43u, 0x6262a6c4u,
    0x9191a839u, 0x9595a431u, 0xe4e437d3u, 0x79798bf2u,
    0xe7e732d5u, 0xc8c8438bu, 0x3737596eu, 0x6d6db7dau,
    0x8d8d8c01u, 0xd5d564b1u, 0x4e4ed29cu, 0xa9a9e049u,
    0x6c6cb4d8u, 0x5656faacu, 0xf4f407f3u, 0xeaea25cfu,
    0x6565afcau, 0x7a7a8ef4u, 0xaeaee947u, 0x08081810u,
    0xbabad56fu, 0x787888f0u, 0x25256f4au, 0x2e2e725cu,
    0x1c1c2438u, 0xa6a6f157u, 0xb4b4c773u, 0xc6c65197u,
    0xe8e823cbu, 0xdddd7ca1u, 0x74749ce8u, 0x1f1f213eu,
    0x4b4bdd96u, 0xbdbddc61u, 0x8b8b860du, 0x8a8a850fu,
    0x707090e0u, 0x3e3e427cu, 0xb5b5c471u, 0x6666aaccu,
    0x4848d890u, 0x03030506u, 0xf6f601f7u, 0x0e0e121cu,
    0x6161a3c2u, 0x35355f6au, 0x5757f9aeu, 0xb9b9d069u,
    0x86869117u, 0xc1c15899u, 0x1d1d273au, 0x9e9eb927u,
    0xe1e138d9u, 0xf8f813ebu, 0x9898b32bu, 0x11113322u,
    0x6969bbd2u, 0xd9d970a9u, 0x8e8e8907u, 0x9494a733u,
    0x9b9bb62du, 0x1e1e223cu, 0x87879215u, 0xe9e920c9u,
    0xcece4987u, 0x5555ffaau, 0x28287850u, 0xdfdf7aa5u,
    0x8c8c8f03u, 0xa1a1f859u, 0x89898009u, 0x0d0d171au,
    0xbfbfda65u, 0xe6e631d7u, 0x4242c684u, 0x6868b8d0u,
    0x4141c382u, 0x9999b029u, 0x2d2d775au, 0x0f0f111eu,
    0xb0b0cb7bu, 0x5454fca8u, 0xbbbbd66du, 0x16163a2cu,
};

#if defined(__riscv) && defined(AES_USE_TT_ASM)
/* Byte-swapped Te tables for the raw little-endian assembly path. */
static const uint32_t aes_te0_r[256] = {
    0xa56363c6u, 0x847c7cf8u, 0x997777eeu, 0x8d7b7bf6u,
    0x0df2f2ffu, 0xbd6b6bd6u, 0xb16f6fdeu, 0x54c5c591u,
    0x50303060u, 0x03010102u, 0xa96767ceu, 0x7d2b2b56u,
    0x19fefee7u, 0x62d7d7b5u, 0xe6abab4du, 0x9a7676ecu,
    0x45caca8fu, 0x9d82821fu, 0x40c9c989u, 0x877d7dfau,
    0x15fafaefu, 0xeb5959b2u, 0xc947478eu, 0x0bf0f0fbu,
    0xecadad41u, 0x67d4d4b3u, 0xfda2a25fu, 0xeaafaf45u,
    0xbf9c9c23u, 0xf7a4a453u, 0x967272e4u, 0x5bc0c09bu,
    0xc2b7b775u, 0x1cfdfde1u, 0xae93933du, 0x6a26264cu,
    0x5a36366cu, 0x413f3f7eu, 0x02f7f7f5u, 0x4fcccc83u,
    0x5c343468u, 0xf4a5a551u, 0x34e5e5d1u, 0x08f1f1f9u,
    0x937171e2u, 0x73d8d8abu, 0x53313162u, 0x3f15152au,
    0x0c040408u, 0x52c7c795u, 0x65232346u, 0x5ec3c39du,
    0x28181830u, 0xa1969637u, 0x0f05050au, 0xb59a9a2fu,
    0x0907070eu, 0x36121224u, 0x9b80801bu, 0x3de2e2dfu,
    0x26ebebcdu, 0x6927274eu, 0xcdb2b27fu, 0x9f7575eau,
    0x1b090912u, 0x9e83831du, 0x742c2c58u, 0x2e1a1a34u,
    0x2d1b1b36u, 0xb26e6edcu, 0xee5a5ab4u, 0xfba0a05bu,
    0xf65252a4u, 0x4d3b3b76u, 0x61d6d6b7u, 0xceb3b37du,
    0x7b292952u, 0x3ee3e3ddu, 0x712f2f5eu, 0x97848413u,
    0xf55353a6u, 0x68d1d1b9u, 0x00000000u, 0x2cededc1u,
    0x60202040u, 0x1ffcfce3u, 0xc8b1b179u, 0xed5b5bb6u,
    0xbe6a6ad4u, 0x46cbcb8du, 0xd9bebe67u, 0x4b393972u,
    0xde4a4a94u, 0xd44c4c98u, 0xe85858b0u, 0x4acfcf85u,
    0x6bd0d0bbu, 0x2aefefc5u, 0xe5aaaa4fu, 0x16fbfbedu,
    0xc5434386u, 0xd74d4d9au, 0x55333366u, 0x94858511u,
    0xcf45458au, 0x10f9f9e9u, 0x06020204u, 0x817f7ffeu,
    0xf05050a0u, 0x443c3c78u, 0xba9f9f25u, 0xe3a8a84bu,
    0xf35151a2u, 0xfea3a35du, 0xc0404080u, 0x8a8f8f05u,
    0xad92923fu, 0xbc9d9d21u, 0x48383870u, 0x04f5f5f1u,
    0xdfbcbc63u, 0xc1b6b677u, 0x75dadaafu, 0x63212142u,
    0x30101020u, 0x1affffe5u, 0x0ef3f3fdu, 0x6dd2d2bfu,
    0x4ccdcd81u, 0x140c0c18u, 0x35131326u, 0x2fececc3u,
    0xe15f5fbeu, 0xa2979735u, 0xcc444488u, 0x3917172eu,
    0x57c4c493u, 0xf2a7a755u, 0x827e7efcu, 0x473d3d7au,
    0xac6464c8u, 0xe75d5dbau, 0x2b191932u, 0x957373e6u,
    0xa06060c0u, 0x98818119u, 0xd14f4f9eu, 0x7fdcdca3u,
    0x66222244u, 0x7e2a2a54u, 0xab90903bu, 0x8388880bu,
    0xca46468cu, 0x29eeeec7u, 0xd3b8b86bu, 0x3c141428u,
    0x79dedea7u, 0xe25e5ebcu, 0x1d0b0b16u, 0x76dbdbadu,
    0x3be0e0dbu, 0x56323264u, 0x4e3a3a74u, 0x1e0a0a14u,
    0xdb494992u, 0x0a06060cu, 0x6c242448u, 0xe45c5cb8u,
    0x5dc2c29fu, 0x6ed3d3bdu, 0xefacac43u, 0xa66262c4u,
    0xa8919139u, 0xa4959531u, 0x37e4e4d3u, 0x8b7979f2u,
    0x32e7e7d5u, 0x43c8c88bu, 0x5937376eu, 0xb76d6ddau,
    0x8c8d8d01u, 0x64d5d5b1u, 0xd24e4e9cu, 0xe0a9a949u,
    0xb46c6cd8u, 0xfa5656acu, 0x07f4f4f3u, 0x25eaeacfu,
    0xaf6565cau, 0x8e7a7af4u, 0xe9aeae47u, 0x18080810u,
    0xd5baba6fu, 0x887878f0u, 0x6f25254au, 0x722e2e5cu,
    0x241c1c38u, 0xf1a6a657u, 0xc7b4b473u, 0x51c6c697u,
    0x23e8e8cbu, 0x7cdddda1u, 0x9c7474e8u, 0x211f1f3eu,
    0xdd4b4b96u, 0xdcbdbd61u, 0x868b8b0du, 0x858a8a0fu,
    0x907070e0u, 0x423e3e7cu, 0xc4b5b571u, 0xaa6666ccu,
    0xd8484890u, 0x05030306u, 0x01f6f6f7u, 0x120e0e1cu,
    0xa36161c2u, 0x5f35356au, 0xf95757aeu, 0xd0b9b969u,
    0x91868617u, 0x58c1c199u, 0x271d1d3au, 0xb99e9e27u,
    0x38e1e1d9u, 0x13f8f8ebu, 0xb398982bu, 0x33111122u,
    0xbb6969d2u, 0x70d9d9a9u, 0x898e8e07u, 0xa7949433u,
    0xb69b9b2du, 0x221e1e3cu, 0x92878715u, 0x20e9e9c9u,
    0x49cece87u, 0xff5555aau, 0x78282850u, 0x7adfdfa5u,
    0x8f8c8c03u, 0xf8a1a159u, 0x80898909u, 0x170d0d1au,
    0xdabfbf65u, 0x31e6e6d7u, 0xc6424284u, 0xb86868d0u,
    0xc3414182u, 0xb0999929u, 0x772d2d5au, 0x110f0f1eu,
    0xcbb0b07bu, 0xfc5454a8u, 0xd6bbbb6du, 0x3a16162cu,
};

static const uint32_t aes_te1_r[256] = {
    0x6363c6a5u, 0x7c7cf884u, 0x7777ee99u, 0x7b7bf68du,
    0xf2f2ff0du, 0x6b6bd6bdu, 0x6f6fdeb1u, 0xc5c59154u,
    0x30306050u, 0x01010203u, 0x6767cea9u, 0x2b2b567du,
    0xfefee719u, 0xd7d7b562u, 0xabab4de6u, 0x7676ec9au,
    0xcaca8f45u, 0x82821f9du, 0xc9c98940u, 0x7d7dfa87u,
    0xfafaef15u, 0x5959b2ebu, 0x47478ec9u, 0xf0f0fb0bu,
    0xadad41ecu, 0xd4d4b367u, 0xa2a25ffdu, 0xafaf45eau,
    0x9c9c23bfu, 0xa4a453f7u, 0x7272e496u, 0xc0c09b5bu,
    0xb7b775c2u, 0xfdfde11cu, 0x93933daeu, 0x26264c6au,
    0x36366c5au, 0x3f3f7e41u, 0xf7f7f502u, 0xcccc834fu,
    0x3434685cu, 0xa5a551f4u, 0xe5e5d134u, 0xf1f1f908u,
    0x7171e293u, 0xd8d8ab73u, 0x31316253u, 0x15152a3fu,
    0x0404080cu, 0xc7c79552u, 0x23234665u, 0xc3c39d5eu,
    0x18183028u, 0x969637a1u, 0x05050a0fu, 0x9a9a2fb5u,
    0x07070e09u, 0x12122436u, 0x80801b9bu, 0xe2e2df3du,
    0xebebcd26u, 0x27274e69u, 0xb2b27fcdu, 0x7575ea9fu,
    0x0909121bu, 0x83831d9eu, 0x2c2c5874u, 0x1a1a342eu,
    0x1b1b362du, 0x6e6edcb2u, 0x5a5ab4eeu, 0xa0a05bfbu,
    0x5252a4f6u, 0x3b3b764du, 0xd6d6b761u, 0xb3b37dceu,
    0x2929527bu, 0xe3e3dd3eu, 0x2f2f5e71u, 0x84841397u,
    0x5353a6f5u, 0xd1d1b968u, 0x00000000u, 0xededc12cu,
    0x20204060u, 0xfcfce31fu, 0xb1b179c8u, 0x5b5bb6edu,
    0x6a6ad4beu, 0xcbcb8d46u, 0xbebe67d9u, 0x3939724bu,
    0x4a4a94deu, 0x4c4c98d4u, 0x5858b0e8u, 0xcfcf854au,
    0xd0d0bb6bu, 0xefefc52au, 0xaaaa4fe5u, 0xfbfbed16u,
    0x434386c5u, 0x4d4d9ad7u, 0x33336655u, 0x85851194u,
    0x45458acfu, 0xf9f9e910u, 0x02020406u, 0x7f7ffe81u,
    0x5050a0f0u, 0x3c3c7844u, 0x9f9f25bau, 0xa8a84be3u,
    0x5151a2f3u, 0xa3a35dfeu, 0x404080c0u, 0x8f8f058au,
    0x92923fadu, 0x9d9d21bcu, 0x38387048u, 0xf5f5f104u,
    0xbcbc63dfu, 0xb6b677c1u, 0xdadaaf75u, 0x21214263u,
    0x10102030u, 0xffffe51au, 0xf3f3fd0eu, 0xd2d2bf6du,
    0xcdcd814cu, 0x0c0c1814u, 0x13132635u, 0xececc32fu,
    0x5f5fbee1u, 0x979735a2u, 0x444488ccu, 0x17172e39u,
    0xc4c49357u, 0xa7a755f2u, 0x7e7efc82u, 0x3d3d7a47u,
    0x6464c8acu, 0x5d5dbae7u, 0x1919322bu, 0x7373e695u,
    0x6060c0a0u, 0x81811998u, 0x4f4f9ed1u, 0xdcdca37fu,
    0x22224466u, 0x2a2a547eu, 0x90903babu, 0x88880b83u,
    0x46468ccau, 0xeeeec729u, 0xb8b86bd3u, 0x1414283cu,
    0xdedea779u, 0x5e5ebce2u, 0x0b0b161du, 0xdbdbad76u,
    0xe0e0db3bu, 0x32326456u, 0x3a3a744eu, 0x0a0a141eu,
    0x494992dbu, 0x06060c0au, 0x2424486cu, 0x5c5cb8e4u,
    0xc2c29f5du, 0xd3d3bd6eu, 0xacac43efu, 0x6262c4a6u,
    0x919139a8u, 0x959531a4u, 0xe4e4d337u, 0x7979f28bu,
    0xe7e7d532u, 0xc8c88b43u, 0x37376e59u, 0x6d6ddab7u,
    0x8d8d018cu, 0xd5d5b164u, 0x4e4e9cd2u, 0xa9a949e0u,
    0x6c6cd8b4u, 0x5656acfau, 0xf4f4f307u, 0xeaeacf25u,
    0x6565caafu, 0x7a7af48eu, 0xaeae47e9u, 0x08081018u,
    0xbaba6fd5u, 0x7878f088u, 0x25254a6fu, 0x2e2e5c72u,
    0x1c1c3824u, 0xa6a657f1u, 0xb4b473c7u, 0xc6c69751u,
    0xe8e8cb23u, 0xdddda17cu, 0x7474e89cu, 0x1f1f3e21u,
    0x4b4b96ddu, 0xbdbd61dcu, 0x8b8b0d86u, 0x8a8a0f85u,
    0x7070e090u, 0x3e3e7c42u, 0xb5b571c4u, 0x6666ccaau,
    0x484890d8u, 0x03030605u, 0xf6f6f701u, 0x0e0e1c12u,
    0x6161c2a3u, 0x35356a5fu, 0x5757aef9u, 0xb9b969d0u,
    0x86861791u, 0xc1c19958u, 0x1d1d3a27u, 0x9e9e27b9u,
    0xe1e1d938u, 0xf8f8eb13u, 0x98982bb3u, 0x11112233u,
    0x6969d2bbu, 0xd9d9a970u, 0x8e8e0789u, 0x949433a7u,
    0x9b9b2db6u, 0x1e1e3c22u, 0x87871592u, 0xe9e9c920u,
    0xcece8749u, 0x5555aaffu, 0x28285078u, 0xdfdfa57au,
    0x8c8c038fu, 0xa1a159f8u, 0x89890980u, 0x0d0d1a17u,
    0xbfbf65dau, 0xe6e6d731u, 0x424284c6u, 0x6868d0b8u,
    0x414182c3u, 0x999929b0u, 0x2d2d5a77u, 0x0f0f1e11u,
    0xb0b07bcbu, 0x5454a8fcu, 0xbbbb6dd6u, 0x16162c3au,
};

static const uint32_t aes_te2_r[256] = {
    0x63c6a563u, 0x7cf8847cu, 0x77ee9977u, 0x7bf68d7bu,
    0xf2ff0df2u, 0x6bd6bd6bu, 0x6fdeb16fu, 0xc59154c5u,
    0x30605030u, 0x01020301u, 0x67cea967u, 0x2b567d2bu,
    0xfee719feu, 0xd7b562d7u, 0xab4de6abu, 0x76ec9a76u,
    0xca8f45cau, 0x821f9d82u, 0xc98940c9u, 0x7dfa877du,
    0xfaef15fau, 0x59b2eb59u, 0x478ec947u, 0xf0fb0bf0u,
    0xad41ecadu, 0xd4b367d4u, 0xa25ffda2u, 0xaf45eaafu,
    0x9c23bf9cu, 0xa453f7a4u, 0x72e49672u, 0xc09b5bc0u,
    0xb775c2b7u, 0xfde11cfdu, 0x933dae93u, 0x264c6a26u,
    0x366c5a36u, 0x3f7e413fu, 0xf7f502f7u, 0xcc834fccu,
    0x34685c34u, 0xa551f4a5u, 0xe5d134e5u, 0xf1f908f1u,
    0x71e29371u, 0xd8ab73d8u, 0x31625331u, 0x152a3f15u,
    0x04080c04u, 0xc79552c7u, 0x23466523u, 0xc39d5ec3u,
    0x18302818u, 0x9637a196u, 0x050a0f05u, 0x9a2fb59au,
    0x070e0907u, 0x12243612u, 0x801b9b80u, 0xe2df3de2u,
    0xebcd26ebu, 0x274e6927u, 0xb27fcdb2u, 0x75ea9f75u,
    0x09121b09u, 0x831d9e83u, 0x2c58742cu, 0x1a342e1au,
    0x1b362d1bu, 0x6edcb26eu, 0x5ab4ee5au, 0xa05bfba0u,
    0x52a4f652u, 0x3b764d3bu, 0xd6b761d6u, 0xb37dceb3u,
    0x29527b29u, 0xe3dd3ee3u, 0x2f5e712fu, 0x84139784u,
    0x53a6f553u, 0xd1b968d1u, 0x00000000u, 0xedc12cedu,
    0x20406020u, 0xfce31ffcu, 0xb179c8b1u, 0x5bb6ed5bu,
    0x6ad4be6au, 0xcb8d46cbu, 0xbe67d9beu, 0x39724b39u,
    0x4a94de4au, 0x4c98d44cu, 0x58b0e858u, 0xcf854acfu,
    0xd0bb6bd0u, 0xefc52aefu, 0xaa4fe5aau, 0xfbed16fbu,
    0x4386c543u, 0x4d9ad74du, 0x33665533u, 0x85119485u,
    0x458acf45u, 0xf9e910f9u, 0x02040602u, 0x7ffe817fu,
    0x50a0f050u, 0x3c78443cu, 0x9f25ba9fu, 0xa84be3a8u,
    0x51a2f351u, 0xa35dfea3u, 0x4080c040u, 0x8f058a8fu,
    0x923fad92u, 0x9d21bc9du, 0x38704838u, 0xf5f104f5u,
    0xbc63dfbcu, 0xb677c1b6u, 0xdaaf75dau, 0x21426321u,
    0x10203010u, 0xffe51affu, 0xf3fd0ef3u, 0xd2bf6dd2u,
    0xcd814ccdu, 0x0c18140cu, 0x13263513u, 0xecc32fecu,
    0x5fbee15fu, 0x9735a297u, 0x4488cc44u, 0x172e3917u,
    0xc49357c4u, 0xa755f2a7u, 0x7efc827eu, 0x3d7a473du,
    0x64c8ac64u, 0x5dbae75du, 0x19322b19u, 0x73e69573u,
    0x60c0a060u, 0x81199881u, 0x4f9ed14fu, 0xdca37fdcu,
    0x22446622u, 0x2a547e2au, 0x903bab90u, 0x880b8388u,
    0x468cca46u, 0xeec729eeu, 0xb86bd3b8u, 0x14283c14u,
    0xdea779deu, 0x5ebce25eu, 0x0b161d0bu, 0xdbad76dbu,
    0xe0db3be0u, 0x32645632u, 0x3a744e3au, 0x0a141e0au,
    0x4992db49u, 0x060c0a06u, 0x24486c24u, 0x5cb8e45cu,
    0xc29f5dc2u, 0xd3bd6ed3u, 0xac43efacu, 0x62c4a662u,
    0x9139a891u, 0x9531a495u, 0xe4d337e4u, 0x79f28b79u,
    0xe7d532e7u, 0xc88b43c8u, 0x376e5937u, 0x6ddab76du,
    0x8d018c8du, 0xd5b164d5u, 0x4e9cd24eu, 0xa949e0a9u,
    0x6cd8b46cu, 0x56acfa56u, 0xf4f307f4u, 0xeacf25eau,
    0x65caaf65u, 0x7af48e7au, 0xae47e9aeu, 0x08101808u,
    0xba6fd5bau, 0x78f08878u, 0x254a6f25u, 0x2e5c722eu,
    0x1c38241cu, 0xa657f1a6u, 0xb473c7b4u, 0xc69751c6u,
    0xe8cb23e8u, 0xdda17cddu, 0x74e89c74u, 0x1f3e211fu,
    0x4b96dd4bu, 0xbd61dcbdu, 0x8b0d868bu, 0x8a0f858au,
    0x70e09070u, 0x3e7c423eu, 0xb571c4b5u, 0x66ccaa66u,
    0x4890d848u, 0x03060503u, 0xf6f701f6u, 0x0e1c120eu,
    0x61c2a361u, 0x356a5f35u, 0x57aef957u, 0xb969d0b9u,
    0x86179186u, 0xc19958c1u, 0x1d3a271du, 0x9e27b99eu,
    0xe1d938e1u, 0xf8eb13f8u, 0x982bb398u, 0x11223311u,
    0x69d2bb69u, 0xd9a970d9u, 0x8e07898eu, 0x9433a794u,
    0x9b2db69bu, 0x1e3c221eu, 0x87159287u, 0xe9c920e9u,
    0xce8749ceu, 0x55aaff55u, 0x28507828u, 0xdfa57adfu,
    0x8c038f8cu, 0xa159f8a1u, 0x89098089u, 0x0d1a170du,
    0xbf65dabfu, 0xe6d731e6u, 0x4284c642u, 0x68d0b868u,
    0x4182c341u, 0x9929b099u, 0x2d5a772du, 0x0f1e110fu,
    0xb07bcbb0u, 0x54a8fc54u, 0xbb6dd6bbu, 0x162c3a16u,
};

static const uint32_t aes_te3_r[256] = {
    0xc6a56363u, 0xf8847c7cu, 0xee997777u, 0xf68d7b7bu,
    0xff0df2f2u, 0xd6bd6b6bu, 0xdeb16f6fu, 0x9154c5c5u,
    0x60503030u, 0x02030101u, 0xcea96767u, 0x567d2b2bu,
    0xe719fefeu, 0xb562d7d7u, 0x4de6ababu, 0xec9a7676u,
    0x8f45cacau, 0x1f9d8282u, 0x8940c9c9u, 0xfa877d7du,
    0xef15fafau, 0xb2eb5959u, 0x8ec94747u, 0xfb0bf0f0u,
    0x41ecadadu, 0xb367d4d4u, 0x5ffda2a2u, 0x45eaafafu,
    0x23bf9c9cu, 0x53f7a4a4u, 0xe4967272u, 0x9b5bc0c0u,
    0x75c2b7b7u, 0xe11cfdfdu, 0x3dae9393u, 0x4c6a2626u,
    0x6c5a3636u, 0x7e413f3fu, 0xf502f7f7u, 0x834fccccu,
    0x685c3434u, 0x51f4a5a5u, 0xd134e5e5u, 0xf908f1f1u,
    0xe2937171u, 0xab73d8d8u, 0x62533131u, 0x2a3f1515u,
    0x080c0404u, 0x9552c7c7u, 0x46652323u, 0x9d5ec3c3u,
    0x30281818u, 0x37a19696u, 0x0a0f0505u, 0x2fb59a9au,
    0x0e090707u, 0x24361212u, 0x1b9b8080u, 0xdf3de2e2u,
    0xcd26ebebu, 0x4e692727u, 0x7fcdb2b2u, 0xea9f7575u,
    0x121b0909u, 0x1d9e8383u, 0x58742c2cu, 0x342e1a1au,
    0x362d1b1bu, 0xdcb26e6eu, 0xb4ee5a5au, 0x5bfba0a0u,
    0xa4f65252u, 0x764d3b3bu, 0xb761d6d6u, 0x7dceb3b3u,
    0x527b2929u, 0xdd3ee3e3u, 0x5e712f2fu, 0x13978484u,
    0xa6f55353u, 0xb968d1d1u, 0x00000000u, 0xc12cededu,
    0x40602020u, 0xe31ffcfcu, 0x79c8b1b1u, 0xb6ed5b5bu,
    0xd4be6a6au, 0x8d46cbcbu, 0x67d9bebeu, 0x724b3939u,
    0x94de4a4au, 0x98d44c4cu, 0xb0e85858u, 0x854acfcfu,
    0xbb6bd0d0u, 0xc52aefefu, 0x4fe5aaaau, 0xed16fbfbu,
    0x86c54343u, 0x9ad74d4du, 0x66553333u, 0x11948585u,
    0x8acf4545u, 0xe910f9f9u, 0x04060202u, 0xfe817f7fu,
    0xa0f05050u, 0x78443c3cu, 0x25ba9f9fu, 0x4be3a8a8u,
    0xa2f35151u, 0x5dfea3a3u, 0x80c04040u, 0x058a8f8fu,
    0x3fad9292u, 0x21bc9d9du, 0x70483838u, 0xf104f5f5u,
    0x63dfbcbcu, 0x77c1b6b6u, 0xaf75dadau, 0x42632121u,
    0x20301010u, 0xe51affffu, 0xfd0ef3f3u, 0xbf6dd2d2u,
    0x814ccdcdu, 0x18140c0cu, 0x26351313u, 0xc32fececu,
    0xbee15f5fu, 0x35a29797u, 0x88cc4444u, 0x2e391717u,
    0x9357c4c4u, 0x55f2a7a7u, 0xfc827e7eu, 0x7a473d3du,
    0xc8ac6464u, 0xbae75d5du, 0x322b1919u, 0xe6957373u,
    0xc0a06060u, 0x19988181u, 0x9ed14f4fu, 0xa37fdcdcu,
    0x44662222u, 0x547e2a2au, 0x3bab9090u, 0x0b838888u,
    0x8cca4646u, 0xc729eeeeu, 0x6bd3b8b8u, 0x283c1414u,
    0xa779dedeu, 0xbce25e5eu, 0x161d0b0bu, 0xad76dbdbu,
    0xdb3be0e0u, 0x64563232u, 0x744e3a3au, 0x141e0a0au,
    0x92db4949u, 0x0c0a0606u, 0x486c2424u, 0xb8e45c5cu,
    0x9f5dc2c2u, 0xbd6ed3d3u, 0x43efacacu, 0xc4a66262u,
    0x39a89191u, 0x31a49595u, 0xd337e4e4u, 0xf28b7979u,
    0xd532e7e7u, 0x8b43c8c8u, 0x6e593737u, 0xdab76d6du,
    0x018c8d8du, 0xb164d5d5u, 0x9cd24e4eu, 0x49e0a9a9u,
    0xd8b46c6cu, 0xacfa5656u, 0xf307f4f4u, 0xcf25eaeau,
    0xcaaf6565u, 0xf48e7a7au, 0x47e9aeaeu, 0x10180808u,
    0x6fd5babau, 0xf0887878u, 0x4a6f2525u, 0x5c722e2eu,
    0x38241c1cu, 0x57f1a6a6u, 0x73c7b4b4u, 0x9751c6c6u,
    0xcb23e8e8u, 0xa17cddddu, 0xe89c7474u, 0x3e211f1fu,
    0x96dd4b4bu, 0x61dcbdbdu, 0x0d868b8bu, 0x0f858a8au,
    0xe0907070u, 0x7c423e3eu, 0x71c4b5b5u, 0xccaa6666u,
    0x90d84848u, 0x06050303u, 0xf701f6f6u, 0x1c120e0eu,
    0xc2a36161u, 0x6a5f3535u, 0xaef95757u, 0x69d0b9b9u,
    0x17918686u, 0x9958c1c1u, 0x3a271d1du, 0x27b99e9eu,
    0xd938e1e1u, 0xeb13f8f8u, 0x2bb39898u, 0x22331111u,
    0xd2bb6969u, 0xa970d9d9u, 0x07898e8eu, 0x33a79494u,
    0x2db69b9bu, 0x3c221e1eu, 0x15928787u, 0xc920e9e9u,
    0x8749ceceu, 0xaaff5555u, 0x50782828u, 0xa57adfdfu,
    0x038f8c8cu, 0x59f8a1a1u, 0x09808989u, 0x1a170d0du,
    0x65dabfbfu, 0xd731e6e6u, 0x84c64242u, 0xd0b86868u,
    0x82c34141u, 0x29b09999u, 0x5a772d2du, 0x1e110f0fu,
    0x7bcbb0b0u, 0xa8fc5454u, 0x6dd6bbbbu, 0x2c3a1616u,
};

#endif
static const uint32_t aes_td0[256] = {
    0x51f4a750u, 0x7e416553u, 0x1a17a4c3u, 0x3a275e96u,
    0x3bab6bcbu, 0x1f9d45f1u, 0xacfa58abu, 0x4be30393u,
    0x2030fa55u, 0xad766df6u, 0x88cc7691u, 0xf5024c25u,
    0x4fe5d7fcu, 0xc52acbd7u, 0x26354480u, 0xb562a38fu,
    0xdeb15a49u, 0x25ba1b67u, 0x45ea0e98u, 0x5dfec0e1u,
    0xc32f7502u, 0x814cf012u, 0x8d4697a3u, 0x6bd3f9c6u,
    0x038f5fe7u, 0x15929c95u, 0xbf6d7aebu, 0x955259dau,
    0xd4be832du, 0x587421d3u, 0x49e06929u, 0x8ec9c844u,
    0x75c2896au, 0xf48e7978u, 0x99583e6bu, 0x27b971ddu,
    0xbee14fb6u, 0xf088ad17u, 0xc920ac66u, 0x7dce3ab4u,
    0x63df4a18u, 0xe51a3182u, 0x97513360u, 0x62537f45u,
    0xb16477e0u, 0xbb6bae84u, 0xfe81a01cu, 0xf9082b94u,
    0x70486858u, 0x8f45fd19u, 0x94de6c87u, 0x527bf8b7u,
    0xab73d323u, 0x724b02e2u, 0xe31f8f57u, 0x6655ab2au,
    0xb2eb2807u, 0x2fb5c203u, 0x86c57b9au, 0xd33708a5u,
    0x302887f2u, 0x23bfa5b2u, 0x02036abau, 0xed16825cu,
    0x8acf1c2bu, 0xa779b492u, 0xf307f2f0u, 0x4e69e2a1u,
    0x65daf4cdu, 0x0605bed5u, 0xd134621fu, 0xc4a6fe8au,
    0x342e539du, 0xa2f355a0u, 0x058ae132u, 0xa4f6eb75u,
    0x0b83ec39u, 0x4060efaau, 0x5e719f06u, 0xbd6e1051u,
    0x3e218af9u, 0x96dd063du, 0xdd3e05aeu, 0x4de6bd46u,
    0x91548db5u, 0x71c45d05u, 0x0406d46fu, 0x605015ffu,
    0x1998fb24u, 0xd6bde997u, 0x894043ccu, 0x67d99e77u,
    0xb0e842bdu, 0x07898b88u, 0xe7195b38u, 0x79c8eedbu,
    0xa17c0a47u, 0x7c420fe9u, 0xf8841ec9u, 0x00000000u,
    0x09808683u, 0x322bed48u, 0x1e1170acu, 0x6c5a724eu,
    0xfd0efffbu, 0x0f853856u, 0x3daed51eu, 0x362d3927u,
    0x0a0fd964u, 0x685ca621u, 0x9b5b54d1u, 0x24362e3au,
    0x0c0a67b1u, 0x9357e70fu, 0xb4ee96d2u, 0x1b9b919eu,
    0x80c0c54fu, 0x61dc20a2u, 0x5a774b69u, 0x1c121a16u,
    0xe293ba0au, 0xc0a02ae5u, 0x3c22e043u, 0x121b171du,
    0x0e090d0bu, 0xf28bc7adu, 0x2db6a8b9u, 0x141ea9c8u,
    0x57f11985u, 0xaf75074cu, 0xee99ddbbu, 0xa37f60fdu,
    0xf701269fu, 0x5c72f5bcu, 0x44663bc5u, 0x5bfb7e34u,
    0x8b432976u, 0xcb23c6dcu, 0xb6edfc68u, 0xb8e4f163u,
    0xd731dccau, 0x42638510u, 0x13972240u, 0x84c61120u,
    0x854a247du, 0xd2bb3df8u, 0xaef93211u, 0xc729a16du,
    0x1d9e2f4bu, 0xdcb230f3u, 0x0d8652ecu, 0x77c1e3d0u,
    0x2bb3166cu, 0xa970b999u, 0x119448fau, 0x47e96422u,
    0xa8fc8cc4u, 0xa0f03f1au, 0x567d2cd8u, 0x223390efu,
    0x87494ec7u, 0xd938d1c1u, 0x8ccaa2feu, 0x98d40b36u,
    0xa6f581cfu, 0xa57ade28u, 0xdab78e26u, 0x3fadbfa4u,
    0x2c3a9de4u, 0x5078920du, 0x6a5fcc9bu, 0x547e4662u,
    0xf68d13c2u, 0x90d8b8e8u, 0x2e39f75eu, 0x82c3aff5u,
    0x9f5d80beu, 0x69d0937cu, 0x6fd52da9u, 0xcf2512b3u,
    0xc8ac993bu, 0x10187da7u, 0xe89c636eu, 0xdb3bbb7bu,
    0xcd267809u, 0x6e5918f4u, 0xec9ab701u, 0x834f9aa8u,
    0xe6956e65u, 0xaaffe67eu, 0x21bccf08u, 0xef15e8e6u,
    0xbae79bd9u, 0x4a6f36ceu, 0xea9f09d4u, 0x29b07cd6u,
    0x31a4b2afu, 0x2a3f2331u, 0xc6a59430u, 0x35a266c0u,
    0x744ebc37u, 0xfc82caa6u, 0xe090d0b0u, 0x33a7d815u,
    0xf104984au, 0x41ecdaf7u, 0x7fcd500eu, 0x1791f62fu,
    0x764dd68du, 0x43efb04du, 0xccaa4d54u, 0xe49604dfu,
    0x9ed1b5e3u, 0x4c6a881bu, 0xc12c1fb8u, 0x4665517fu,
    0x9d5eea04u, 0x018c355du, 0xfa877473u, 0xfb0b412eu,
    0xb3671d5au, 0x92dbd252u, 0xe9105633u, 0x6dd64713u,
    0x9ad7618cu, 0x37a10c7au, 0x59f8148eu, 0xeb133c89u,
    0xcea927eeu, 0xb761c935u, 0xe11ce5edu, 0x7a47b13cu,
    0x9cd2df59u, 0x55f2733fu, 0x1814ce79u, 0x73c737bfu,
    0x53f7cdeau, 0x5ffdaa5bu, 0xdf3d6f14u, 0x7844db86u,
    0xcaaff381u, 0xb968c43eu, 0x3824342cu, 0xc2a3405fu,
    0x161dc372u, 0xbce2250cu, 0x283c498bu, 0xff0d9541u,
    0x39a80171u, 0x080cb3deu, 0xd8b4e49cu, 0x6456c190u,
    0x7bcb8461u, 0xd532b670u, 0x486c5c74u, 0xd0b85742u,
};

static const uint32_t aes_td1[256] = {
    0x5051f4a7u, 0x537e4165u, 0xc31a17a4u, 0x963a275eu,
    0xcb3bab6bu, 0xf11f9d45u, 0xabacfa58u, 0x934be303u,
    0x552030fau, 0xf6ad766du, 0x9188cc76u, 0x25f5024cu,
    0xfc4fe5d7u, 0xd7c52acbu, 0x80263544u, 0x8fb562a3u,
    0x49deb15au, 0x6725ba1bu, 0x9845ea0eu, 0xe15dfec0u,
    0x02c32f75u, 0x12814cf0u, 0xa38d4697u, 0xc66bd3f9u,
    0xe7038f5fu, 0x9515929cu, 0xebbf6d7au, 0xda955259u,
    0x2dd4be83u, 0xd3587421u, 0x2949e069u, 0x448ec9c8u,
    0x6a75c289u, 0x78f48e79u, 0x6b99583eu, 0xdd27b971u,
    0xb6bee14fu, 0x17f088adu, 0x66c920acu, 0xb47dce3au,
    0x1863df4au, 0x82e51a31u, 0x60975133u, 0x4562537fu,
    0xe0b16477u, 0x84bb6baeu, 0x1cfe81a0u, 0x94f9082bu,
    0x58704868u, 0x198f45fdu, 0x8794de6cu, 0xb7527bf8u,
    0x23ab73d3u, 0xe2724b02u, 0x57e31f8fu, 0x2a6655abu,
    0x07b2eb28u, 0x032fb5c2u, 0x9a86c57bu, 0xa5d33708u,
    0xf2302887u, 0xb223bfa5u, 0xba02036au, 0x5ced1682u,
    0x2b8acf1cu, 0x92a779b4u, 0xf0f307f2u, 0xa14e69e2u,
    0xcd65daf4u, 0xd50605beu, 0x1fd13462u, 0x8ac4a6feu,
    0x9d342e53u, 0xa0a2f355u, 0x32058ae1u, 0x75a4f6ebu,
    0x390b83ecu, 0xaa4060efu, 0x065e719fu, 0x51bd6e10u,
    0xf93e218au, 0x3d96dd06u, 0xaedd3e05u, 0x464de6bdu,
    0xb591548du, 0x0571c45du, 0x6f0406d4u, 0xff605015u,
    0x241998fbu, 0x97d6bde9u, 0xcc894043u, 0x7767d99eu,
    0xbdb0e842u, 0x8807898bu, 0x38e7195bu, 0xdb79c8eeu,
    0x47a17c0au, 0xe97c420fu, 0xc9f8841eu, 0x00000000u,
    0x83098086u, 0x48322bedu, 0xac1e1170u, 0x4e6c5a72u,
    0xfbfd0effu, 0x560f8538u, 0x1e3daed5u, 0x27362d39u,
    0x640a0fd9u, 0x21685ca6u, 0xd19b5b54u, 0x3a24362eu,
    0xb10c0a67u, 0x0f9357e7u, 0xd2b4ee96u, 0x9e1b9b91u,
    0x4f80c0c5u, 0xa261dc20u, 0x695a774bu, 0x161c121au,
    0x0ae293bau, 0xe5c0a02au, 0x433c22e0u, 0x1d121b17u,
    0x0b0e090du, 0xadf28bc7u, 0xb92db6a8u, 0xc8141ea9u,
    0x8557f119u, 0x4caf7507u, 0xbbee99ddu, 0xfda37f60u,
    0x9ff70126u, 0xbc5c72f5u, 0xc544663bu, 0x345bfb7eu,
    0x768b4329u, 0xdccb23c6u, 0x68b6edfcu, 0x63b8e4f1u,
    0xcad731dcu, 0x10426385u, 0x40139722u, 0x2084c611u,
    0x7d854a24u, 0xf8d2bb3du, 0x11aef932u, 0x6dc729a1u,
    0x4b1d9e2fu, 0xf3dcb230u, 0xec0d8652u, 0xd077c1e3u,
    0x6c2bb316u, 0x99a970b9u, 0xfa119448u, 0x2247e964u,
    0xc4a8fc8cu, 0x1aa0f03fu, 0xd8567d2cu, 0xef223390u,
    0xc787494eu, 0xc1d938d1u, 0xfe8ccaa2u, 0x3698d40bu,
    0xcfa6f581u, 0x28a57adeu, 0x26dab78eu, 0xa43fadbfu,
    0xe42c3a9du, 0x0d507892u, 0x9b6a5fccu, 0x62547e46u,
    0xc2f68d13u, 0xe890d8b8u, 0x5e2e39f7u, 0xf582c3afu,
    0xbe9f5d80u, 0x7c69d093u, 0xa96fd52du, 0xb3cf2512u,
    0x3bc8ac99u, 0xa710187du, 0x6ee89c63u, 0x7bdb3bbbu,
    0x09cd2678u, 0xf46e5918u, 0x01ec9ab7u, 0xa8834f9au,
    0x65e6956eu, 0x7eaaffe6u, 0x0821bccfu, 0xe6ef15e8u,
    0xd9bae79bu, 0xce4a6f36u, 0xd4ea9f09u, 0xd629b07cu,
    0xaf31a4b2u, 0x312a3f23u, 0x30c6a594u, 0xc035a266u,
    0x37744ebcu, 0xa6fc82cau, 0xb0e090d0u, 0x1533a7d8u,
    0x4af10498u, 0xf741ecdau, 0x0e7fcd50u, 0x2f1791f6u,
    0x8d764dd6u, 0x4d43efb0u, 0x54ccaa4du, 0xdfe49604u,
    0xe39ed1b5u, 0x1b4c6a88u, 0xb8c12c1fu, 0x7f466551u,
    0x049d5eeau, 0x5d018c35u, 0x73fa8774u, 0x2efb0b41u,
    0x5ab3671du, 0x5292dbd2u, 0x33e91056u, 0x136dd647u,
    0x8c9ad761u, 0x7a37a10cu, 0x8e59f814u, 0x89eb133cu,
    0xeecea927u, 0x35b761c9u, 0xede11ce5u, 0x3c7a47b1u,
    0x599cd2dfu, 0x3f55f273u, 0x791814ceu, 0xbf73c737u,
    0xea53f7cdu, 0x5b5ffdaau, 0x14df3d6fu, 0x867844dbu,
    0x81caaff3u, 0x3eb968c4u, 0x2c382434u, 0x5fc2a340u,
    0x72161dc3u, 0x0cbce225u, 0x8b283c49u, 0x41ff0d95u,
    0x7139a801u, 0xde080cb3u, 0x9cd8b4e4u, 0x906456c1u,
    0x617bcb84u, 0x70d532b6u, 0x74486c5cu, 0x42d0b857u,
};

static const uint32_t aes_td2[256] = {
    0xa75051f4u, 0x65537e41u, 0xa4c31a17u, 0x5e963a27u,
    0x6bcb3babu, 0x45f11f9du, 0x58abacfau, 0x03934be3u,
    0xfa552030u, 0x6df6ad76u, 0x769188ccu, 0x4c25f502u,
    0xd7fc4fe5u, 0xcbd7c52au, 0x44802635u, 0xa38fb562u,
    0x5a49deb1u, 0x1b6725bau, 0x0e9845eau, 0xc0e15dfeu,
    0x7502c32fu, 0xf012814cu, 0x97a38d46u, 0xf9c66bd3u,
    0x5fe7038fu, 0x9c951592u, 0x7aebbf6du, 0x59da9552u,
    0x832dd4beu, 0x21d35874u, 0x692949e0u, 0xc8448ec9u,
    0x896a75c2u, 0x7978f48eu, 0x3e6b9958u, 0x71dd27b9u,
    0x4fb6bee1u, 0xad17f088u, 0xac66c920u, 0x3ab47dceu,
    0x4a1863dfu, 0x3182e51au, 0x33609751u, 0x7f456253u,
    0x77e0b164u, 0xae84bb6bu, 0xa01cfe81u, 0x2b94f908u,
    0x68587048u, 0xfd198f45u, 0x6c8794deu, 0xf8b7527bu,
    0xd323ab73u, 0x02e2724bu, 0x8f57e31fu, 0xab2a6655u,
    0x2807b2ebu, 0xc2032fb5u, 0x7b9a86c5u, 0x08a5d337u,
    0x87f23028u, 0xa5b223bfu, 0x6aba0203u, 0x825ced16u,
    0x1c2b8acfu, 0xb492a779u, 0xf2f0f307u, 0xe2a14e69u,
    0xf4cd65dau, 0xbed50605u, 0x621fd134u, 0xfe8ac4a6u,
    0x539d342eu, 0x55a0a2f3u, 0xe132058au, 0xeb75a4f6u,
    0xec390b83u, 0xefaa4060u, 0x9f065e71u, 0x1051bd6eu,
    0x8af93e21u, 0x063d96ddu, 0x05aedd3eu, 0xbd464de6u,
    0x8db59154u, 0x5d0571c4u, 0xd46f0406u, 0x15ff6050u,
    0xfb241998u, 0xe997d6bdu, 0x43cc8940u, 0x9e7767d9u,
    0x42bdb0e8u, 0x8b880789u, 0x5b38e719u, 0xeedb79c8u,
    0x0a47a17cu, 0x0fe97c42u, 0x1ec9f884u, 0x00000000u,
    0x86830980u, 0xed48322bu, 0x70ac1e11u, 0x724e6c5au,
    0xfffbfd0eu, 0x38560f85u, 0xd51e3daeu, 0x3927362du,
    0xd9640a0fu, 0xa621685cu, 0x54d19b5bu, 0x2e3a2436u,
    0x67b10c0au, 0xe70f9357u, 0x96d2b4eeu, 0x919e1b9bu,
    0xc54f80c0u, 0x20a261dcu, 0x4b695a77u, 0x1a161c12u,
    0xba0ae293u, 0x2ae5c0a0u, 0xe0433c22u, 0x171d121bu,
    0x0d0b0e09u, 0xc7adf28bu, 0xa8b92db6u, 0xa9c8141eu,
    0x198557f1u, 0x074caf75u, 0xddbbee99u, 0x60fda37fu,
    0x269ff701u, 0xf5bc5c72u, 0x3bc54466u, 0x7e345bfbu,
    0x29768b43u, 0xc6dccb23u, 0xfc68b6edu, 0xf163b8e4u,
    0xdccad731u, 0x85104263u, 0x22401397u, 0x112084c6u,
    0x247d854au, 0x3df8d2bbu, 0x3211aef9u, 0xa16dc729u,
    0x2f4b1d9eu, 0x30f3dcb2u, 0x52ec0d86u, 0xe3d077c1u,
    0x166c2bb3u, 0xb999a970u, 0x48fa1194u, 0x642247e9u,
    0x8cc4a8fcu, 0x3f1aa0f0u, 0x2cd8567du, 0x90ef2233u,
    0x4ec78749u, 0xd1c1d938u, 0xa2fe8ccau, 0x0b3698d4u,
    0x81cfa6f5u, 0xde28a57au, 0x8e26dab7u, 0xbfa43fadu,
    0x9de42c3au, 0x920d5078u, 0xcc9b6a5fu, 0x4662547eu,
    0x13c2f68du, 0xb8e890d8u, 0xf75e2e39u, 0xaff582c3u,
    0x80be9f5du, 0x937c69d0u, 0x2da96fd5u, 0x12b3cf25u,
    0x993bc8acu, 0x7da71018u, 0x636ee89cu, 0xbb7bdb3bu,
    0x7809cd26u, 0x18f46e59u, 0xb701ec9au, 0x9aa8834fu,
    0x6e65e695u, 0xe67eaaffu, 0xcf0821bcu, 0xe8e6ef15u,
    0x9bd9bae7u, 0x36ce4a6fu, 0x09d4ea9fu, 0x7cd629b0u,
    0xb2af31a4u, 0x23312a3fu, 0x9430c6a5u, 0x66c035a2u,
    0xbc37744eu, 0xcaa6fc82u, 0xd0b0e090u, 0xd81533a7u,
    0x984af104u, 0xdaf741ecu, 0x500e7fcdu, 0xf62f1791u,
    0xd68d764du, 0xb04d43efu, 0x4d54ccaau, 0x04dfe496u,
    0xb5e39ed1u, 0x881b4c6au, 0x1fb8c12cu, 0x517f4665u,
    0xea049d5eu, 0x355d018cu, 0x7473fa87u, 0x412efb0bu,
    0x1d5ab367u, 0xd25292dbu, 0x5633e910u, 0x47136dd6u,
    0x618c9ad7u, 0x0c7a37a1u, 0x148e59f8u, 0x3c89eb13u,
    0x27eecea9u, 0xc935b761u, 0xe5ede11cu, 0xb13c7a47u,
    0xdf599cd2u, 0x733f55f2u, 0xce791814u, 0x37bf73c7u,
    0xcdea53f7u, 0xaa5b5ffdu, 0x6f14df3du, 0xdb867844u,
    0xf381caafu, 0xc43eb968u, 0x342c3824u, 0x405fc2a3u,
    0xc372161du, 0x250cbce2u, 0x498b283cu, 0x9541ff0du,
    0x017139a8u, 0xb3de080cu, 0xe49cd8b4u, 0xc1906456u,
    0x84617bcbu, 0xb670d532u, 0x5c74486cu, 0x5742d0b8u,
};

static const uint32_t aes_td3[256] = {
    0xf4a75051u, 0x4165537eu, 0x17a4c31au, 0x275e963au,
    0xab6bcb3bu, 0x9d45f11fu, 0xfa58abacu, 0xe303934bu,
    0x30fa5520u, 0x766df6adu, 0xcc769188u, 0x024c25f5u,
    0xe5d7fc4fu, 0x2acbd7c5u, 0x35448026u, 0x62a38fb5u,
    0xb15a49deu, 0xba1b6725u, 0xea0e9845u, 0xfec0e15du,
    0x2f7502c3u, 0x4cf01281u, 0x4697a38du, 0xd3f9c66bu,
    0x8f5fe703u, 0x929c9515u, 0x6d7aebbfu, 0x5259da95u,
    0xbe832dd4u, 0x7421d358u, 0xe0692949u, 0xc9c8448eu,
    0xc2896a75u, 0x8e7978f4u, 0x583e6b99u, 0xb971dd27u,
    0xe14fb6beu, 0x88ad17f0u, 0x20ac66c9u, 0xce3ab47du,
    0xdf4a1863u, 0x1a3182e5u, 0x51336097u, 0x537f4562u,
    0x6477e0b1u, 0x6bae84bbu, 0x81a01cfeu, 0x082b94f9u,
    0x48685870u, 0x45fd198fu, 0xde6c8794u, 0x7bf8b752u,
    0x73d323abu, 0x4b02e272u, 0x1f8f57e3u, 0x55ab2a66u,
    0xeb2807b2u, 0xb5c2032fu, 0xc57b9a86u, 0x3708a5d3u,
    0x2887f230u, 0xbfa5b223u, 0x036aba02u, 0x16825cedu,
    0xcf1c2b8au, 0x79b492a7u, 0x07f2f0f3u, 0x69e2a14eu,
    0xdaf4cd65u, 0x05bed506u, 0x34621fd1u, 0xa6fe8ac4u,
    0x2e539d34u, 0xf355a0a2u, 0x8ae13205u, 0xf6eb75a4u,
    0x83ec390bu, 0x60efaa40u, 0x719f065eu, 0x6e1051bdu,
    0x218af93eu, 0xdd063d96u, 0x3e05aeddu, 0xe6bd464du,
    0x548db591u, 0xc45d0571u, 0x06d46f04u, 0x5015ff60u,
    0x98fb2419u, 0xbde997d6u, 0x4043cc89u, 0xd99e7767u,
    0xe842bdb0u, 0x898b8807u, 0x195b38e7u, 0xc8eedb79u,
    0x7c0a47a1u, 0x420fe97cu, 0x841ec9f8u, 0x00000000u,
    0x80868309u, 0x2bed4832u, 0x1170ac1eu, 0x5a724e6cu,
    0x0efffbfdu, 0x8538560fu, 0xaed51e3du, 0x2d392736u,
    0x0fd9640au, 0x5ca62168u, 0x5b54d19bu, 0x362e3a24u,
    0x0a67b10cu, 0x57e70f93u, 0xee96d2b4u, 0x9b919e1bu,
    0xc0c54f80u, 0xdc20a261u, 0x774b695au, 0x121a161cu,
    0x93ba0ae2u, 0xa02ae5c0u, 0x22e0433cu, 0x1b171d12u,
    0x090d0b0eu, 0x8bc7adf2u, 0xb6a8b92du, 0x1ea9c814u,
    0xf1198557u, 0x75074cafu, 0x99ddbbeeu, 0x7f60fda3u,
    0x01269ff7u, 0x72f5bc5cu, 0x663bc544u, 0xfb7e345bu,
    0x4329768bu, 0x23c6dccbu, 0xedfc68b6u, 0xe4f163b8u,
    0x31dccad7u, 0x63851042u, 0x97224013u, 0xc6112084u,
    0x4a247d85u, 0xbb3df8d2u, 0xf93211aeu, 0x29a16dc7u,
    0x9e2f4b1du, 0xb230f3dcu, 0x8652ec0du, 0xc1e3d077u,
    0xb3166c2bu, 0x70b999a9u, 0x9448fa11u, 0xe9642247u,
    0xfc8cc4a8u, 0xf03f1aa0u, 0x7d2cd856u, 0x3390ef22u,
    0x494ec787u, 0x38d1c1d9u, 0xcaa2fe8cu, 0xd40b3698u,
    0xf581cfa6u, 0x7ade28a5u, 0xb78e26dau, 0xadbfa43fu,
    0x3a9de42cu, 0x78920d50u, 0x5fcc9b6au, 0x7e466254u,
    0x8d13c2f6u, 0xd8b8e890u, 0x39f75e2eu, 0xc3aff582u,
    0x5d80be9fu, 0xd0937c69u, 0xd52da96fu, 0x2512b3cfu,
    0xac993bc8u, 0x187da710u, 0x9c636ee8u, 0x3bbb7bdbu,
    0x267809cdu, 0x5918f46eu, 0x9ab701ecu, 0x4f9aa883u,
    0x956e65e6u, 0xffe67eaau, 0xbccf0821u, 0x15e8e6efu,
    0xe79bd9bau, 0x6f36ce4au, 0x9f09d4eau, 0xb07cd629u,
    0xa4b2af31u, 0x3f23312au, 0xa59430c6u, 0xa266c035u,
    0x4ebc3774u, 0x82caa6fcu, 0x90d0b0e0u, 0xa7d81533u,
    0x04984af1u, 0xecdaf741u, 0xcd500e7fu, 0x91f62f17u,
    0x4dd68d76u, 0xefb04d43u, 0xaa4d54ccu, 0x9604dfe4u,
    0xd1b5e39eu, 0x6a881b4cu, 0x2c1fb8c1u, 0x65517f46u,
    0x5eea049du, 0x8c355d01u, 0x877473fau, 0x0b412efbu,
    0x671d5ab3u, 0xdbd25292u, 0x105633e9u, 0xd647136du,
    0xd7618c9au, 0xa10c7a37u, 0xf8148e59u, 0x133c89ebu,
    0xa927eeceu, 0x61c935b7u, 0x1ce5ede1u, 0x47b13c7au,
    0xd2df599cu, 0xf2733f55u, 0x14ce7918u, 0xc737bf73u,
    0xf7cdea53u, 0xfdaa5b5fu, 0x3d6f14dfu, 0x44db8678u,
    0xaff381cau, 0x68c43eb9u, 0x24342c38u, 0xa3405fc2u,
    0x1dc37216u, 0xe2250cbcu, 0x3c498b28u, 0x0d9541ffu,
    0xa8017139u, 0x0cb3de08u, 0xb4e49cd8u, 0x56c19064u,
    0xcb84617bu, 0x32b670d5u, 0x6c5c7448u, 0xb85742d0u,
};

/*
 * Bitsliced byte state.  Eight bit planes are packed four per 64-bit word:
 * plane k lives at bits [16*k, 16*k+15] of p[k/4], and bit i of a plane is
 * bit i of byte i.  Packing four 16-bit lanes per register keeps the
 * fixed GF(2^8) operations cheap on RV64 without any table access.
 */
#define BS_MASK16 0xffffu

static uint64_t bs_plane(const uint64_t p[2], int k)
{
    return (p[k >> 2] >> (16 * (k & 3))) & BS_MASK16;
}

/* Transpose an 8x8 bit matrix held in a little-endian 64-bit word. */
static uint64_t bs_transpose8(uint64_t x)
{
    uint64_t t;

    t = (x ^ (x >> 7)) & 0x00aa00aa00aa00aaULL;
    x ^= t ^ (t << 7);
    t = (x ^ (x >> 14)) & 0x0000cccc0000ccccULL;
    x ^= t ^ (t << 14);
    t = (x ^ (x >> 28)) & 0x00000000f0f0f0f0ULL;
    x ^= t ^ (t << 28);
    return x;
}

static void bs_from_bytes(uint64_t p[2], const uint8_t *s, int n)
{
    uint8_t tmp[16];
    uint64_t x0;
    uint64_t x1;
    uint64_t t0;
    uint64_t t1;
    int i;
    int k;

    for (i = 0; i < n; i++) {
        tmp[i] = s[i];
    }
    for (; i < 16; i++) {
        tmp[i] = 0;
    }
    memcpy(&x0, tmp, sizeof(x0));
    memcpy(&x1, tmp + 8, sizeof(x1));
    t0 = bs_transpose8(x0);
    t1 = bs_transpose8(x1);

    p[0] = 0;
    p[1] = 0;
    for (k = 0; k < 8; k++) {
        uint64_t plane = ((t0 >> (8 * k)) & 0xffu) |
                         (((t1 >> (8 * k)) & 0xffu) << 8);
        if (k < 4) {
            p[0] |= plane << (16 * k);
        } else {
            p[1] |= plane << (16 * (k - 4));
        }
    }
}

static void bs_to_bytes(uint8_t *s, const uint64_t p[2], int n)
{
    uint8_t tmp[16];
    uint64_t t0 = 0;
    uint64_t t1 = 0;
    uint64_t x0;
    uint64_t x1;
    int i;
    int k;

    for (k = 0; k < 8; k++) {
        uint64_t plane = bs_plane(p, k);

        t0 |= (plane & 0xffu) << (8 * k);
        t1 |= ((plane >> 8) & 0xffu) << (8 * k);
    }
    x0 = bs_transpose8(t0);
    x1 = bs_transpose8(t1);
    memcpy(tmp, &x0, sizeof(x0));
    memcpy(tmp + 8, &x1, sizeof(x1));
    for (i = 0; i < n; i++) {
        s[i] = tmp[i];
    }
}

/* Schoolbook GF(2^8) multiplication across sixteen independent lanes. */
static void bs_mul(uint64_t r[2], const uint64_t a[2], const uint64_t b[2])
{
    uint64_t t0 = a[0];
    uint64_t t1 = a[1];
    uint64_t acc0 = 0;
    uint64_t acc1 = 0;

#define BS_MUL_STEP(i)                                                     \
    do {                                                                   \
        uint64_t bi = bs_plane(b, (i));                                    \
        uint64_t m = bi | (bi << 16) | (bi << 32) | (bi << 48);            \
        uint64_t hi;                                                       \
        uint64_t n0;                                                       \
        uint64_t n1;                                                       \
                                                                           \
        acc0 ^= t0 & m;                                                    \
        acc1 ^= t1 & m;                                                    \
        hi = t1 >> 48;                                                     \
        n0 = t0 << 16;                                                     \
        n1 = (t1 << 16) | (t0 >> 48);                                      \
        n0 ^= hi | (hi << 16) | (hi << 48);                                \
        n1 ^= hi;                                                          \
        t0 = n0;                                                           \
        t1 = n1;                                                           \
    } while (0)

    BS_MUL_STEP(0);
    BS_MUL_STEP(1);
    BS_MUL_STEP(2);
    BS_MUL_STEP(3);
    BS_MUL_STEP(4);
    BS_MUL_STEP(5);
    BS_MUL_STEP(6);
    BS_MUL_STEP(7);

#undef BS_MUL_STEP

    r[0] = acc0;
    r[1] = acc1;
}

/* Squaring is linear in GF(2^8): use the fixed basis map instead of a mul. */
static void bs_square(uint64_t r[2], const uint64_t a[2])
{
    uint64_t a0 = bs_plane(a, 0);
    uint64_t a1 = bs_plane(a, 1);
    uint64_t a2 = bs_plane(a, 2);
    uint64_t a3 = bs_plane(a, 3);
    uint64_t a4 = bs_plane(a, 4);
    uint64_t a5 = bs_plane(a, 5);
    uint64_t a6 = bs_plane(a, 6);
    uint64_t a7 = bs_plane(a, 7);
    uint64_t r0;
    uint64_t r1;

    r0 = (a0 ^ a4 ^ a6) |
         ((a4 ^ a6 ^ a7) << 16) |
         ((a1 ^ a5) << 32) |
         ((a4 ^ a5 ^ a6 ^ a7) << 48);
    r1 = (a2 ^ a4 ^ a7) |
         ((a5 ^ a6) << 16) |
         ((a3 ^ a5) << 32) |
         ((a6 ^ a7) << 48);
    r[0] = r0;
    r[1] = r1;
}

/* x^254 = x^-1, computed by a fixed addition chain of length 11:
 * 1,2,4,8,16,32,64,80,84,168,252,254. */
static void bs_inverse(uint64_t r[2], const uint64_t a[2])
{
    uint64_t x2[2], x4[2], x8[2], x16[2], x32[2], x64[2];
    uint64_t x80[2], x84[2], x168[2], x252[2];

    bs_square(x2, a);
    bs_square(x4, x2);
    bs_square(x8, x4);
    bs_square(x16, x8);
    bs_square(x32, x16);
    bs_square(x64, x32);

    bs_mul(x80, x64, x16);
    bs_mul(x84, x80, x4);
    bs_square(x168, x84);
    bs_mul(x252, x168, x84);
    bs_mul(r, x252, x2);
}

/* AES affine transform: x ^ rotl(x,1) ^ rotl(x,2) ^ rotl(x,3) ^
 * rotl(x,4) ^ 0x63.  Rotations are plane permutations. */
static void bs_affine(uint64_t r[2], const uint64_t a[2])
{
    uint64_t a0 = bs_plane(a, 0);
    uint64_t a1 = bs_plane(a, 1);
    uint64_t a2 = bs_plane(a, 2);
    uint64_t a3 = bs_plane(a, 3);
    uint64_t a4 = bs_plane(a, 4);
    uint64_t a5 = bs_plane(a, 5);
    uint64_t a6 = bs_plane(a, 6);
    uint64_t a7 = bs_plane(a, 7);

    r[0] = (a0 ^ a7 ^ a6 ^ a5 ^ a4 ^ BS_MASK16) |
           ((a1 ^ a0 ^ a7 ^ a6 ^ a5 ^ BS_MASK16) << 16) |
           ((a2 ^ a1 ^ a0 ^ a7 ^ a6) << 32) |
           ((a3 ^ a2 ^ a1 ^ a0 ^ a7) << 48);
    r[1] = (a4 ^ a3 ^ a2 ^ a1 ^ a0) |
           ((a5 ^ a4 ^ a3 ^ a2 ^ a1 ^ BS_MASK16) << 16) |
           ((a6 ^ a5 ^ a4 ^ a3 ^ a2 ^ BS_MASK16) << 32) |
           ((a7 ^ a6 ^ a5 ^ a4 ^ a3) << 48);
}

/* Inverse affine transform, followed by the GF inverse in inv_sbox. */
static void bs_inverse_affine(uint64_t r[2], const uint64_t a[2])
{
    uint64_t a0 = bs_plane(a, 0);
    uint64_t a1 = bs_plane(a, 1);
    uint64_t a2 = bs_plane(a, 2);
    uint64_t a3 = bs_plane(a, 3);
    uint64_t a4 = bs_plane(a, 4);
    uint64_t a5 = bs_plane(a, 5);
    uint64_t a6 = bs_plane(a, 6);
    uint64_t a7 = bs_plane(a, 7);

    r[0] = (a7 ^ a5 ^ a2 ^ BS_MASK16) |
           ((a0 ^ a6 ^ a3) << 16) |
           ((a1 ^ a7 ^ a4 ^ BS_MASK16) << 32) |
           ((a2 ^ a0 ^ a5) << 48);
    r[1] = (a3 ^ a1 ^ a6) |
           ((a4 ^ a2 ^ a7) << 16) |
           ((a5 ^ a3 ^ a0) << 32) |
           ((a6 ^ a4 ^ a1) << 48);
}

static void bs_sbox(uint64_t p[2])
{
    uint64_t t[2];

    bs_inverse(t, p);
    bs_affine(p, t);
}

static void bs_inv_sbox(uint64_t p[2])
{
    uint64_t t[2];

    bs_inverse_affine(t, p);
    bs_inverse(p, t);
}

#if defined(__riscv) && defined(AES_USE_ZVKNED)
static int g_use_asm = 1;
#else
static int g_use_asm = 0;
#endif
static int g_use_fast;

#if defined(__riscv) && defined(AES_USE_RVV)
static int g_use_rvv = 0;

int aes_set_rvv(int use_rvv)
{
    if (use_rvv && !aes_rvv_available()) {
        return -1;
    }
    g_use_rvv = use_rvv ? 1 : 0;
    if (g_use_rvv) {
        g_use_asm = 0;
        g_use_fast = 0;
    }
    return 0;
}

int aes_get_rvv(void)
{
    return g_use_rvv;
}
#else
int aes_set_rvv(int use_rvv)
{
    return use_rvv ? -1 : 0;
}

int aes_get_rvv(void)
{
    return 0;
}
#endif

int aes_set_impl(int use_asm)
{
    if (use_asm == 2) {
        g_use_fast = 1;
        g_use_asm = 0;
#if defined(__riscv) && defined(AES_USE_RVV)
        g_use_rvv = 0;
#endif
        return 0;
    }
    g_use_fast = 0;
#if defined(__riscv) && defined(AES_USE_ZVKNED)
    if (use_asm && !g_use_asm) {
        g_use_asm = 1;
#if defined(__riscv) && defined(AES_USE_RVV)
        g_use_rvv = 0;
#endif
    } else if (!use_asm) {
        g_use_asm = 0;
    }
    return 0;
#else
    (void)use_asm;
    g_use_asm = 0;
    return use_asm ? -1 : 0;
#endif
}

int aes_get_impl(void)
{
    if (g_use_fast) {
        return 3;
    }
    if (g_use_asm) {
        return 1;
    }
#if defined(__riscv) && defined(AES_USE_RVV)
    if (g_use_rvv) {
        return 2;
    }
#endif
    return 0;
}

static uint32_t load_be32(const uint8_t *p)
{
    uint32_t v;

    memcpy(&v, p, sizeof(v));
    return __builtin_bswap32(v);
}

static void store_be32(uint8_t *p, uint32_t v)
{
    v = __builtin_bswap32(v);
    memcpy(p, &v, sizeof(v));
}

static uint32_t bswap32(uint32_t v)
{
    return ((v & 0x000000ffu) << 24) | ((v & 0x0000ff00u) << 8) |
           ((v & 0x00ff0000u) >> 8) | ((v & 0xff000000u) >> 24);
}

static uint8_t gmul(uint8_t a, uint8_t b);
static uint32_t inv_mix_word(uint32_t w);
static void aes_encrypt_block_tt(uint8_t out[16], const uint8_t in[16],
                                 const aes_key *ctx);
static void aes_decrypt_block_tt(uint8_t out[16], const uint8_t in[16],
                                 const aes_key *ctx);
#if defined(__riscv) && defined(AES_USE_TT_ASM)
void aes_tt_rounds_rv64(uint32_t s[4], const uint32_t *rk, int nrounds,
                        const uint32_t *te0, const uint32_t *te1,
                        const uint32_t *te2, const uint32_t *te3);
void aes_encrypt_block_tt_rv64(uint8_t out[16], const uint8_t in[16],
                               const uint32_t *rk, int rounds,
                               const uint32_t *te0, const uint32_t *te1,
                               const uint32_t *te2, const uint32_t *te3,
                               const uint8_t *sbox);
#endif
static void aes_encrypt4_tt(uint8_t out[64], const uint8_t in[64],
                            const aes_key *ctx);
static void aes_decrypt4_tt(uint8_t out[64], const uint8_t in[64],
                            const aes_key *ctx);

static uint32_t sub_word(uint32_t w)
{
    uint8_t b[4];
    uint64_t p[2];

    b[0] = (uint8_t)(w >> 24);
    b[1] = (uint8_t)(w >> 16);
    b[2] = (uint8_t)(w >> 8);
    b[3] = (uint8_t)w;
    bs_from_bytes(p, b, 4);
    bs_sbox(p);
    bs_to_bytes(b, p, 4);
    return ((uint32_t)b[0] << 24) | ((uint32_t)b[1] << 16) |
           ((uint32_t)b[2] << 8) | (uint32_t)b[3];
}

static uint32_t inv_mix_word(uint32_t w)
{
    uint8_t a0 = (uint8_t)(w >> 24);
    uint8_t a1 = (uint8_t)(w >> 16);
    uint8_t a2 = (uint8_t)(w >> 8);
    uint8_t a3 = (uint8_t)w;

    return ((uint32_t)(gmul(a0, 14) ^ gmul(a1, 11) ^
                       gmul(a2, 13) ^ gmul(a3, 9)) << 24) |
           ((uint32_t)(gmul(a0, 9) ^ gmul(a1, 14) ^
                       gmul(a2, 11) ^ gmul(a3, 13)) << 16) |
           ((uint32_t)(gmul(a0, 13) ^ gmul(a1, 9) ^
                       gmul(a2, 14) ^ gmul(a3, 11)) << 8) |
           (uint32_t)(gmul(a0, 11) ^ gmul(a1, 13) ^
                      gmul(a2, 9) ^ gmul(a3, 14));
}

static uint32_t rot_word(uint32_t w)
{
    return (w << 8) | (w >> 24);
}

int aes_set_encrypt_key(const uint8_t *key, int bits, aes_key *ctx)
{
    int nk;
    int nr;
    int total;
    int i;

    if (bits == 128) {
        nk = 4;
    } else if (bits == 192) {
        nk = 6;
    } else if (bits == 256) {
        nk = 8;
    } else {
        return -1;
    }
    nr = nk + 6;
    total = 4 * (nr + 1);
    if (total > AES_MAX_ROUND_KEYS) {
        return -1;
    }

    for (i = 0; i < nk; i++) {
        ctx->rk[i] = load_be32(key + 4 * i);
    }
    for (i = nk; i < total; i++) {
        uint32_t temp = ctx->rk[i - 1];

        if (i % nk == 0) {
            temp = sub_word(rot_word(temp)) ^ ((uint32_t)aes_rcon[i / nk] << 24);
        } else if (nk > 6 && i % nk == 4) {
            temp = sub_word(temp);
        }
        ctx->rk[i] = ctx->rk[i - nk] ^ temp;
    }
    for (i = 0; i < total; i++) {
        ctx->rk_z[i] = bswap32(ctx->rk[i]);
        ctx->rk_b[4 * i + 0] = (uint8_t)(ctx->rk[i] >> 24);
        ctx->rk_b[4 * i + 1] = (uint8_t)(ctx->rk[i] >> 16);
        ctx->rk_b[4 * i + 2] = (uint8_t)(ctx->rk[i] >> 8);
        ctx->rk_b[4 * i + 3] = (uint8_t)ctx->rk[i];
    }
    for (i = 0; i < 4; i++) {
        ctx->rk_d[i] = ctx->rk[4 * nr + i];
        ctx->rk_d[4 * nr + i] = ctx->rk[i];
    }
    for (i = 1; i < nr; i++) {
        int c;
        for (c = 0; c < 4; c++) {
            ctx->rk_d[4 * (nr - i) + c] = inv_mix_word(ctx->rk[4 * i + c]);
        }
    }
    ctx->rounds = nr;
    return 0;
}

static uint8_t xtime(uint8_t x)
{
    uint8_t hi = (uint8_t)(0u - (uint8_t)(x >> 7));
    return (uint8_t)((x << 1) ^ (hi & 0x1b));
}

static void add_round_key(uint8_t s[16], const uint32_t *rk, int round)
{
    int c;
    int r;

    for (c = 0; c < 4; c++) {
        uint32_t w = rk[4 * round + c];
        for (r = 0; r < 4; r++) {
            s[r + 4 * c] ^= (uint8_t)(w >> (24 - 8 * r));
        }
    }
}

static void sub_bytes(uint8_t s[16])
{
    uint64_t p[2];

    bs_from_bytes(p, s, 16);
    bs_sbox(p);
    bs_to_bytes(s, p, 16);
}

static void shift_rows(uint8_t s[16])
{
    uint8_t t[16];
    int r;
    int c;

    for (r = 0; r < 4; r++) {
        for (c = 0; c < 4; c++) {
            t[r + 4 * c] = s[r + 4 * ((c + r) & 3)];
        }
    }
    memcpy(s, t, sizeof(t));
}

static void mix_columns(uint8_t s[16])
{
    int c;

    for (c = 0; c < 4; c++) {
        uint8_t *p = s + 4 * c;
        uint8_t a0 = p[0], a1 = p[1], a2 = p[2], a3 = p[3];
        p[0] = (uint8_t)(xtime(a0) ^ (xtime(a1) ^ a1) ^ a2 ^ a3);
        p[1] = (uint8_t)(a0 ^ xtime(a1) ^ (xtime(a2) ^ a2) ^ a3);
        p[2] = (uint8_t)(a0 ^ a1 ^ xtime(a2) ^ (xtime(a3) ^ a3));
        p[3] = (uint8_t)((xtime(a0) ^ a0) ^ a1 ^ a2 ^ xtime(a3));
    }
}

void aes_encrypt_block_c(uint8_t out[16], const uint8_t in[16],
                         const aes_key *ctx)
{
    uint8_t s[16];
    int round;

    memcpy(s, in, sizeof(s));
    add_round_key(s, ctx->rk, 0);
    for (round = 1; round < ctx->rounds; round++) {
        sub_bytes(s);
        shift_rows(s);
        mix_columns(s);
        add_round_key(s, ctx->rk, round);
    }
    sub_bytes(s);
    shift_rows(s);
    add_round_key(s, ctx->rk, ctx->rounds);
    memcpy(out, s, sizeof(s));
}

void aes_encrypt_block(uint8_t out[16], const uint8_t in[16],
                       const aes_key *ctx)
{
    if (g_use_fast) {
        aes_encrypt_block_tt(out, in, ctx);
        return;
    }
#if defined(__riscv) && defined(AES_USE_ZVKNED)
    if (g_use_asm) {
        aes_encrypt_block_zvkned(out, in, ctx->rk_z, ctx->rounds);
        return;
    }
#endif
#if defined(__riscv) && defined(AES_USE_RVV)
    if (g_use_rvv) {
        aes_encrypt_block_rvv(out, in, ctx->rk_b, ctx->rounds);
        return;
    }
#endif
    aes_encrypt_block_c(out, in, ctx);
}

void aes_ecb_encrypt(uint8_t *out, const uint8_t *in, size_t len,
                     const aes_key *ctx)
{
    size_t off;

#if defined(__riscv) && defined(AES_USE_ZVKNED)
    if (g_use_asm && ctx->rounds == 10) {
        aes_ecb_encrypt_zvkned_128(out, in, len, ctx->rk_z);
        return;
    }
    if (g_use_asm && ctx->rounds == 12) {
        aes_ecb_encrypt_zvkned_192(out, in, len, ctx->rk_z);
        return;
    }
    if (g_use_asm && ctx->rounds == 14) {
        aes_ecb_encrypt_zvkned_256(out, in, len, ctx->rk_z);
        return;
    }
#endif
    if (g_use_fast) {
        for (off = 0; off + 16 <= len; off += 16) {
            aes_encrypt_block_tt(out + off, in + off, ctx);
        }
        return;
    }
    for (off = 0; off + 16 <= len; off += 16) {
        aes_encrypt_block(out + off, in + off, ctx);
    }
}

void aes_cbc_encrypt(uint8_t *out, const uint8_t *in, size_t len,
                     uint8_t iv[16], const aes_key *ctx)
{
    uint8_t block[16];
    size_t off;
    int i;

#if defined(__riscv) && defined(AES_USE_ZVKNED)
    if (g_use_asm && ctx->rounds == 10) {
        aes_cbc_encrypt_zvkned_128(out, in, len, iv, ctx->rk_z);
        return;
    }
    if (g_use_asm && ctx->rounds == 12) {
        aes_cbc_encrypt_zvkned_192(out, in, len, iv, ctx->rk_z);
        return;
    }
    if (g_use_asm && ctx->rounds == 14) {
        aes_cbc_encrypt_zvkned_256(out, in, len, iv, ctx->rk_z);
        return;
    }
#endif
    if (g_use_fast) {
        for (off = 0; off + 16 <= len; off += 16) {
            uint64_t a;
            uint64_t b;
            uint64_t c;
            uint64_t d;

            memcpy(&a, in + off, 8);
            memcpy(&b, in + off + 8, 8);
            memcpy(&c, iv, 8);
            memcpy(&d, iv + 8, 8);
            a ^= c;
            b ^= d;
            memcpy(block, &a, 8);
            memcpy(block + 8, &b, 8);
            aes_encrypt_block_tt(out + off, block, ctx);
            memcpy(iv, out + off, 16);
        }
        return;
    }
    for (off = 0; off + 16 <= len; off += 16) {
        for (i = 0; i < 16; i++) {
            block[i] = (uint8_t)(in[off + i] ^ iv[i]);
        }
        aes_encrypt_block(out + off, block, ctx);
        memcpy(iv, out + off, 16);
    }
}

static uint8_t gmul(uint8_t a, uint8_t b)
{
    uint8_t result = 0;
    int i;

    for (i = 0; i < 8; i++) {
        uint8_t mask = (uint8_t)(0u - (uint8_t)(b & 1u));
        result ^= a & mask;
        a = xtime(a);
        b >>= 1;
    }
    return result;
}

static void aes_encrypt_block_tt(uint8_t out[16], const uint8_t in[16],
                                 const aes_key *ctx)
{
    uint32_t s0;
    uint32_t s1;
    uint32_t s2;
    uint32_t s3;

#if defined(__riscv) && defined(AES_USE_TT_ASM)
    if ((((uintptr_t)in | (uintptr_t)out) & 3u) == 0) {
        aes_encrypt_block_tt_rv64(out, in, ctx->rk_z, ctx->rounds,
                                  aes_te0_r, aes_te1_r, aes_te2_r, aes_te3_r,
                                  aes_sbox_tt);
        return;
    }
#endif
    {
        uint32_t t0, t1, t2, t3;

        s0 = load_be32(in) ^ ctx->rk[0];
        s1 = load_be32(in + 4) ^ ctx->rk[1];
        s2 = load_be32(in + 8) ^ ctx->rk[2];
        s3 = load_be32(in + 12) ^ ctx->rk[3];

#define TT_ROUND(a0, a1, a2, a3, b0, b1, b2, b3, k)                        \
        do {                                                                \
            (b0) = aes_te0[(a0) >> 24] ^                                    \
                   aes_te1[((a1) >> 16) & 0xff] ^                           \
                   aes_te2[((a2) >> 8) & 0xff] ^                            \
                   aes_te3[(a3) & 0xff] ^ (k)[0];                           \
            (b1) = aes_te0[(a1) >> 24] ^                                    \
                   aes_te1[((a2) >> 16) & 0xff] ^                           \
                   aes_te2[((a3) >> 8) & 0xff] ^                            \
                   aes_te3[(a0) & 0xff] ^ (k)[1];                           \
            (b2) = aes_te0[(a2) >> 24] ^                                    \
                   aes_te1[((a3) >> 16) & 0xff] ^                           \
                   aes_te2[((a0) >> 8) & 0xff] ^                            \
                   aes_te3[(a1) & 0xff] ^ (k)[2];                           \
            (b3) = aes_te0[(a3) >> 24] ^                                    \
                   aes_te1[((a0) >> 16) & 0xff] ^                           \
                   aes_te2[((a1) >> 8) & 0xff] ^                            \
                   aes_te3[(a2) & 0xff] ^ (k)[3];                           \
        } while (0)

        TT_ROUND(s0, s1, s2, s3, t0, t1, t2, t3, ctx->rk + 4);
        TT_ROUND(t0, t1, t2, t3, s0, s1, s2, s3, ctx->rk + 8);
        TT_ROUND(s0, s1, s2, s3, t0, t1, t2, t3, ctx->rk + 12);
        TT_ROUND(t0, t1, t2, t3, s0, s1, s2, s3, ctx->rk + 16);
        TT_ROUND(s0, s1, s2, s3, t0, t1, t2, t3, ctx->rk + 20);
        TT_ROUND(t0, t1, t2, t3, s0, s1, s2, s3, ctx->rk + 24);
        TT_ROUND(s0, s1, s2, s3, t0, t1, t2, t3, ctx->rk + 28);
        TT_ROUND(t0, t1, t2, t3, s0, s1, s2, s3, ctx->rk + 32);
        TT_ROUND(s0, s1, s2, s3, t0, t1, t2, t3, ctx->rk + 36);
        if (ctx->rounds > 10) {
            TT_ROUND(t0, t1, t2, t3, s0, s1, s2, s3, ctx->rk + 40);
            TT_ROUND(s0, s1, s2, s3, t0, t1, t2, t3, ctx->rk + 44);
            if (ctx->rounds > 12) {
                TT_ROUND(t0, t1, t2, t3, s0, s1, s2, s3, ctx->rk + 48);
                TT_ROUND(s0, s1, s2, s3, t0, t1, t2, t3, ctx->rk + 52);
            }
        }
#undef TT_ROUND

        /* The unrolled chain always finishes in t0..t3 (odd round count). */
        s0 = t0;
        s1 = t1;
        s2 = t2;
        s3 = t3;
    }
    {
        const uint32_t *rk = ctx->rk + 4 * ctx->rounds;
        uint32_t w0 = ((uint32_t)aes_sbox_tt[s0 >> 24] << 24) |
                      ((uint32_t)aes_sbox_tt[(s1 >> 16) & 0xff] << 16) |
                      ((uint32_t)aes_sbox_tt[(s2 >> 8) & 0xff] << 8) |
                      (uint32_t)aes_sbox_tt[s3 & 0xff];
        uint32_t w1 = ((uint32_t)aes_sbox_tt[s1 >> 24] << 24) |
                      ((uint32_t)aes_sbox_tt[(s2 >> 16) & 0xff] << 16) |
                      ((uint32_t)aes_sbox_tt[(s3 >> 8) & 0xff] << 8) |
                      (uint32_t)aes_sbox_tt[s0 & 0xff];
        uint32_t w2 = ((uint32_t)aes_sbox_tt[s2 >> 24] << 24) |
                      ((uint32_t)aes_sbox_tt[(s3 >> 16) & 0xff] << 16) |
                      ((uint32_t)aes_sbox_tt[(s0 >> 8) & 0xff] << 8) |
                      (uint32_t)aes_sbox_tt[s1 & 0xff];
        uint32_t w3 = ((uint32_t)aes_sbox_tt[s3 >> 24] << 24) |
                      ((uint32_t)aes_sbox_tt[(s0 >> 16) & 0xff] << 16) |
                      ((uint32_t)aes_sbox_tt[(s1 >> 8) & 0xff] << 8) |
                      (uint32_t)aes_sbox_tt[s2 & 0xff];

        store_be32(out, w0 ^ rk[0]);
        store_be32(out + 4, w1 ^ rk[1]);
        store_be32(out + 8, w2 ^ rk[2]);
        store_be32(out + 12, w3 ^ rk[3]);
    }
}

static void aes_decrypt_block_tt(uint8_t out[16], const uint8_t in[16],
                                 const aes_key *ctx)
{
    uint32_t s0 = load_be32(in) ^ ctx->rk_d[0];
    uint32_t s1 = load_be32(in + 4) ^ ctx->rk_d[1];
    uint32_t s2 = load_be32(in + 8) ^ ctx->rk_d[2];
    uint32_t s3 = load_be32(in + 12) ^ ctx->rk_d[3];
    uint32_t t0, t1, t2, t3;

#define TT_DROUND(a0, a1, a2, a3, b0, b1, b2, b3, k)                            do {                                                                            (b0) = aes_td0[(a0) >> 24] ^                                                       aes_td1[((a3) >> 16) & 0xff] ^                                              aes_td2[((a2) >> 8) & 0xff] ^                                               aes_td3[(a1) & 0xff] ^ (k)[0];                                       (b1) = aes_td0[(a1) >> 24] ^                                                       aes_td1[((a0) >> 16) & 0xff] ^                                              aes_td2[((a3) >> 8) & 0xff] ^                                               aes_td3[(a2) & 0xff] ^ (k)[1];                                       (b2) = aes_td0[(a2) >> 24] ^                                                       aes_td1[((a1) >> 16) & 0xff] ^                                              aes_td2[((a0) >> 8) & 0xff] ^                                               aes_td3[(a3) & 0xff] ^ (k)[2];                                       (b3) = aes_td0[(a3) >> 24] ^                                                       aes_td1[((a2) >> 16) & 0xff] ^                                              aes_td2[((a1) >> 8) & 0xff] ^                                               aes_td3[(a0) & 0xff] ^ (k)[3];                                   } while (0)

    TT_DROUND(s0, s1, s2, s3, t0, t1, t2, t3, ctx->rk_d + 4);
    TT_DROUND(t0, t1, t2, t3, s0, s1, s2, s3, ctx->rk_d + 8);
    TT_DROUND(s0, s1, s2, s3, t0, t1, t2, t3, ctx->rk_d + 12);
    TT_DROUND(t0, t1, t2, t3, s0, s1, s2, s3, ctx->rk_d + 16);
    TT_DROUND(s0, s1, s2, s3, t0, t1, t2, t3, ctx->rk_d + 20);
    TT_DROUND(t0, t1, t2, t3, s0, s1, s2, s3, ctx->rk_d + 24);
    TT_DROUND(s0, s1, s2, s3, t0, t1, t2, t3, ctx->rk_d + 28);
    TT_DROUND(t0, t1, t2, t3, s0, s1, s2, s3, ctx->rk_d + 32);
    TT_DROUND(s0, s1, s2, s3, t0, t1, t2, t3, ctx->rk_d + 36);
    if (ctx->rounds > 10) {
        TT_DROUND(t0, t1, t2, t3, s0, s1, s2, s3, ctx->rk_d + 40);
        TT_DROUND(s0, s1, s2, s3, t0, t1, t2, t3, ctx->rk_d + 44);
        if (ctx->rounds > 12) {
            TT_DROUND(t0, t1, t2, t3, s0, s1, s2, s3, ctx->rk_d + 48);
            TT_DROUND(s0, s1, s2, s3, t0, t1, t2, t3, ctx->rk_d + 52);
        }
    }
#undef TT_DROUND
    s0 = t0;
    s1 = t1;
    s2 = t2;
    s3 = t3;
    {
        const uint32_t *rk = ctx->rk_d + 4 * ctx->rounds;
        uint32_t w0 = ((uint32_t)aes_inv_sbox_tt[s0 >> 24] << 24) |
                      ((uint32_t)aes_inv_sbox_tt[(s3 >> 16) & 0xff] << 16) |
                      ((uint32_t)aes_inv_sbox_tt[(s2 >> 8) & 0xff] << 8) |
                      (uint32_t)aes_inv_sbox_tt[s1 & 0xff];
        uint32_t w1 = ((uint32_t)aes_inv_sbox_tt[s1 >> 24] << 24) |
                      ((uint32_t)aes_inv_sbox_tt[(s0 >> 16) & 0xff] << 16) |
                      ((uint32_t)aes_inv_sbox_tt[(s3 >> 8) & 0xff] << 8) |
                      (uint32_t)aes_inv_sbox_tt[s2 & 0xff];
        uint32_t w2 = ((uint32_t)aes_inv_sbox_tt[s2 >> 24] << 24) |
                      ((uint32_t)aes_inv_sbox_tt[(s1 >> 16) & 0xff] << 16) |
                      ((uint32_t)aes_inv_sbox_tt[(s0 >> 8) & 0xff] << 8) |
                      (uint32_t)aes_inv_sbox_tt[s3 & 0xff];
        uint32_t w3 = ((uint32_t)aes_inv_sbox_tt[s3 >> 24] << 24) |
                      ((uint32_t)aes_inv_sbox_tt[(s2 >> 16) & 0xff] << 16) |
                      ((uint32_t)aes_inv_sbox_tt[(s1 >> 8) & 0xff] << 8) |
                      (uint32_t)aes_inv_sbox_tt[s0 & 0xff];

        store_be32(out, w0 ^ rk[0]);
        store_be32(out + 4, w1 ^ rk[1]);
        store_be32(out + 8, w2 ^ rk[2]);
        store_be32(out + 12, w3 ^ rk[3]);
    }
}

static void inv_sub_bytes(uint8_t s[16])
{
    uint64_t p[2];

    bs_from_bytes(p, s, 16);
    bs_inv_sbox(p);
    bs_to_bytes(s, p, 16);
}

static void inv_shift_rows(uint8_t s[16])
{
    uint8_t t[16];
    int r;
    int c;

    for (r = 0; r < 4; r++) {
        for (c = 0; c < 4; c++) {
            t[r + 4 * c] = s[r + 4 * ((c - r + 4) & 3)];
        }
    }
    memcpy(s, t, sizeof(t));
}

/* Four-block T-table kernels.  Interleaving independent blocks hides the
 * table-load latency on in-order cores. */
#define TT_ENC_STEP(x0, x1, x2, x3, k0, k1, k2, k3)                        \
    do {                                                                    \
        uint32_t u0 = aes_te0[(x0) >> 24] ^                                \
                      aes_te1[((x1) >> 16) & 0xff] ^                        \
                      aes_te2[((x2) >> 8) & 0xff] ^                         \
                      aes_te3[(x3) & 0xff] ^ (k0);                          \
        uint32_t u1 = aes_te0[(x1) >> 24] ^                                \
                      aes_te1[((x2) >> 16) & 0xff] ^                        \
                      aes_te2[((x3) >> 8) & 0xff] ^                         \
                      aes_te3[(x0) & 0xff] ^ (k1);                          \
        uint32_t u2 = aes_te0[(x2) >> 24] ^                                \
                      aes_te1[((x3) >> 16) & 0xff] ^                        \
                      aes_te2[((x0) >> 8) & 0xff] ^                         \
                      aes_te3[(x1) & 0xff] ^ (k2);                          \
        uint32_t u3 = aes_te0[(x3) >> 24] ^                                \
                      aes_te1[((x0) >> 16) & 0xff] ^                        \
                      aes_te2[((x1) >> 8) & 0xff] ^                         \
                      aes_te3[(x2) & 0xff] ^ (k3);                          \
        (x0) = u0;                                                         \
        (x1) = u1;                                                         \
        (x2) = u2;                                                         \
        (x3) = u3;                                                         \
    } while (0)

#define TT_ENC_FINAL(x0, x1, x2, x3, off, k0, k1, k2, k3)                  \
    do {                                                                    \
        uint32_t w0 = ((uint32_t)aes_sbox_tt[(x0) >> 24] << 24) |          \
                      ((uint32_t)aes_sbox_tt[((x1) >> 16) & 0xff] << 16) |  \
                      ((uint32_t)aes_sbox_tt[((x2) >> 8) & 0xff] << 8) |    \
                      (uint32_t)aes_sbox_tt[(x3) & 0xff];                   \
        uint32_t w1 = ((uint32_t)aes_sbox_tt[(x1) >> 24] << 24) |          \
                      ((uint32_t)aes_sbox_tt[((x2) >> 16) & 0xff] << 16) |  \
                      ((uint32_t)aes_sbox_tt[((x3) >> 8) & 0xff] << 8) |    \
                      (uint32_t)aes_sbox_tt[(x0) & 0xff];                   \
        uint32_t w2 = ((uint32_t)aes_sbox_tt[(x2) >> 24] << 24) |          \
                      ((uint32_t)aes_sbox_tt[((x3) >> 16) & 0xff] << 16) |  \
                      ((uint32_t)aes_sbox_tt[((x0) >> 8) & 0xff] << 8) |    \
                      (uint32_t)aes_sbox_tt[(x1) & 0xff];                   \
        uint32_t w3 = ((uint32_t)aes_sbox_tt[(x3) >> 24] << 24) |          \
                      ((uint32_t)aes_sbox_tt[((x0) >> 16) & 0xff] << 16) |  \
                      ((uint32_t)aes_sbox_tt[((x1) >> 8) & 0xff] << 8) |    \
                      (uint32_t)aes_sbox_tt[(x2) & 0xff];                   \
        store_be32((out) + (off) + 0, w0 ^ (k0));                          \
        store_be32((out) + (off) + 4, w1 ^ (k1));                          \
        store_be32((out) + (off) + 8, w2 ^ (k2));                          \
        store_be32((out) + (off) + 12, w3 ^ (k3));                         \
    } while (0)

#define TT_DEC_STEP(x0, x1, x2, x3, k0, k1, k2, k3)                        \
    do {                                                                    \
        uint32_t u0 = aes_td0[(x0) >> 24] ^                                \
                      aes_td1[((x3) >> 16) & 0xff] ^                        \
                      aes_td2[((x2) >> 8) & 0xff] ^                         \
                      aes_td3[(x1) & 0xff] ^ (k0);                          \
        uint32_t u1 = aes_td0[(x1) >> 24] ^                                \
                      aes_td1[((x0) >> 16) & 0xff] ^                        \
                      aes_td2[((x3) >> 8) & 0xff] ^                         \
                      aes_td3[(x2) & 0xff] ^ (k1);                          \
        uint32_t u2 = aes_td0[(x2) >> 24] ^                                \
                      aes_td1[((x1) >> 16) & 0xff] ^                        \
                      aes_td2[((x0) >> 8) & 0xff] ^                         \
                      aes_td3[(x3) & 0xff] ^ (k2);                          \
        uint32_t u3 = aes_td0[(x3) >> 24] ^                                \
                      aes_td1[((x2) >> 16) & 0xff] ^                        \
                      aes_td2[((x1) >> 8) & 0xff] ^                         \
                      aes_td3[(x0) & 0xff] ^ (k3);                          \
        (x0) = u0;                                                         \
        (x1) = u1;                                                         \
        (x2) = u2;                                                         \
        (x3) = u3;                                                         \
    } while (0)

#define TT_DEC_FINAL(x0, x1, x2, x3, off, k0, k1, k2, k3)                  \
    do {                                                                    \
        uint32_t w0 = ((uint32_t)aes_inv_sbox_tt[(x0) >> 24] << 24) |      \
                      ((uint32_t)aes_inv_sbox_tt[((x3) >> 16) & 0xff] << 16) | \
                      ((uint32_t)aes_inv_sbox_tt[((x2) >> 8) & 0xff] << 8) | \
                      (uint32_t)aes_inv_sbox_tt[(x1) & 0xff];               \
        uint32_t w1 = ((uint32_t)aes_inv_sbox_tt[(x1) >> 24] << 24) |      \
                      ((uint32_t)aes_inv_sbox_tt[((x0) >> 16) & 0xff] << 16) | \
                      ((uint32_t)aes_inv_sbox_tt[((x3) >> 8) & 0xff] << 8) | \
                      (uint32_t)aes_inv_sbox_tt[(x2) & 0xff];               \
        uint32_t w2 = ((uint32_t)aes_inv_sbox_tt[(x2) >> 24] << 24) |      \
                      ((uint32_t)aes_inv_sbox_tt[((x1) >> 16) & 0xff] << 16) | \
                      ((uint32_t)aes_inv_sbox_tt[((x0) >> 8) & 0xff] << 8) | \
                      (uint32_t)aes_inv_sbox_tt[(x3) & 0xff];               \
        uint32_t w3 = ((uint32_t)aes_inv_sbox_tt[(x3) >> 24] << 24) |      \
                      ((uint32_t)aes_inv_sbox_tt[((x2) >> 16) & 0xff] << 16) | \
                      ((uint32_t)aes_inv_sbox_tt[((x1) >> 8) & 0xff] << 8) | \
                      (uint32_t)aes_inv_sbox_tt[(x0) & 0xff];               \
        store_be32((out) + (off) + 0, w0 ^ (k0));                          \
        store_be32((out) + (off) + 4, w1 ^ (k1));                          \
        store_be32((out) + (off) + 8, w2 ^ (k2));                          \
        store_be32((out) + (off) + 12, w3 ^ (k3));                         \
    } while (0)

static void aes_encrypt4_tt(uint8_t out[64], const uint8_t in[64],
                            const aes_key *ctx)
{
    uint32_t a0 = load_be32(in + 0) ^ ctx->rk[0];
    uint32_t a1 = load_be32(in + 4) ^ ctx->rk[1];
    uint32_t a2 = load_be32(in + 8) ^ ctx->rk[2];
    uint32_t a3 = load_be32(in + 12) ^ ctx->rk[3];
    uint32_t b0 = load_be32(in + 16) ^ ctx->rk[0];
    uint32_t b1 = load_be32(in + 20) ^ ctx->rk[1];
    uint32_t b2 = load_be32(in + 24) ^ ctx->rk[2];
    uint32_t b3 = load_be32(in + 28) ^ ctx->rk[3];
    uint32_t c0 = load_be32(in + 32) ^ ctx->rk[0];
    uint32_t c1 = load_be32(in + 36) ^ ctx->rk[1];
    uint32_t c2 = load_be32(in + 40) ^ ctx->rk[2];
    uint32_t c3 = load_be32(in + 44) ^ ctx->rk[3];
    uint32_t d0 = load_be32(in + 48) ^ ctx->rk[0];
    uint32_t d1 = load_be32(in + 52) ^ ctx->rk[1];
    uint32_t d2 = load_be32(in + 56) ^ ctx->rk[2];
    uint32_t d3 = load_be32(in + 60) ^ ctx->rk[3];
    int round;

    for (round = 1; round < ctx->rounds; round++) {
        const uint32_t *rk = ctx->rk + 4 * round;

        TT_ENC_STEP(a0, a1, a2, a3, rk[0], rk[1], rk[2], rk[3]);
        TT_ENC_STEP(b0, b1, b2, b3, rk[0], rk[1], rk[2], rk[3]);
        TT_ENC_STEP(c0, c1, c2, c3, rk[0], rk[1], rk[2], rk[3]);
        TT_ENC_STEP(d0, d1, d2, d3, rk[0], rk[1], rk[2], rk[3]);
    }
    {
        const uint32_t *rk = ctx->rk + 4 * ctx->rounds;

        TT_ENC_FINAL(a0, a1, a2, a3, 0, rk[0], rk[1], rk[2], rk[3]);
        TT_ENC_FINAL(b0, b1, b2, b3, 16, rk[0], rk[1], rk[2], rk[3]);
        TT_ENC_FINAL(c0, c1, c2, c3, 32, rk[0], rk[1], rk[2], rk[3]);
        TT_ENC_FINAL(d0, d1, d2, d3, 48, rk[0], rk[1], rk[2], rk[3]);
    }
}

static void aes_decrypt4_tt(uint8_t out[64], const uint8_t in[64],
                            const aes_key *ctx)
{
    uint32_t a0 = load_be32(in + 0) ^ ctx->rk_d[0];
    uint32_t a1 = load_be32(in + 4) ^ ctx->rk_d[1];
    uint32_t a2 = load_be32(in + 8) ^ ctx->rk_d[2];
    uint32_t a3 = load_be32(in + 12) ^ ctx->rk_d[3];
    uint32_t b0 = load_be32(in + 16) ^ ctx->rk_d[0];
    uint32_t b1 = load_be32(in + 20) ^ ctx->rk_d[1];
    uint32_t b2 = load_be32(in + 24) ^ ctx->rk_d[2];
    uint32_t b3 = load_be32(in + 28) ^ ctx->rk_d[3];
    uint32_t c0 = load_be32(in + 32) ^ ctx->rk_d[0];
    uint32_t c1 = load_be32(in + 36) ^ ctx->rk_d[1];
    uint32_t c2 = load_be32(in + 40) ^ ctx->rk_d[2];
    uint32_t c3 = load_be32(in + 44) ^ ctx->rk_d[3];
    uint32_t d0 = load_be32(in + 48) ^ ctx->rk_d[0];
    uint32_t d1 = load_be32(in + 52) ^ ctx->rk_d[1];
    uint32_t d2 = load_be32(in + 56) ^ ctx->rk_d[2];
    uint32_t d3 = load_be32(in + 60) ^ ctx->rk_d[3];
    int round;

    for (round = 1; round < ctx->rounds; round++) {
        const uint32_t *rk = ctx->rk_d + 4 * round;

        TT_DEC_STEP(a0, a1, a2, a3, rk[0], rk[1], rk[2], rk[3]);
        TT_DEC_STEP(b0, b1, b2, b3, rk[0], rk[1], rk[2], rk[3]);
        TT_DEC_STEP(c0, c1, c2, c3, rk[0], rk[1], rk[2], rk[3]);
        TT_DEC_STEP(d0, d1, d2, d3, rk[0], rk[1], rk[2], rk[3]);
    }
    {
        const uint32_t *rk = ctx->rk_d + 4 * ctx->rounds;

        TT_DEC_FINAL(a0, a1, a2, a3, 0, rk[0], rk[1], rk[2], rk[3]);
        TT_DEC_FINAL(b0, b1, b2, b3, 16, rk[0], rk[1], rk[2], rk[3]);
        TT_DEC_FINAL(c0, c1, c2, c3, 32, rk[0], rk[1], rk[2], rk[3]);
        TT_DEC_FINAL(d0, d1, d2, d3, 48, rk[0], rk[1], rk[2], rk[3]);
    }
}

static void inv_mix_columns(uint8_t s[16])
{
    static const uint8_t m[4] = { 14, 11, 13, 9 };
    int c;

    for (c = 0; c < 4; c++) {
        uint8_t *p = s + 4 * c;
        uint8_t a0 = p[0], a1 = p[1], a2 = p[2], a3 = p[3];
        p[0] = (uint8_t)(gmul(a0, m[0]) ^ gmul(a1, m[1]) ^ gmul(a2, m[2]) ^ gmul(a3, m[3]));
        p[1] = (uint8_t)(gmul(a0, m[3]) ^ gmul(a1, m[0]) ^ gmul(a2, m[1]) ^ gmul(a3, m[2]));
        p[2] = (uint8_t)(gmul(a0, m[2]) ^ gmul(a1, m[3]) ^ gmul(a2, m[0]) ^ gmul(a3, m[1]));
        p[3] = (uint8_t)(gmul(a0, m[1]) ^ gmul(a1, m[2]) ^ gmul(a2, m[3]) ^ gmul(a3, m[0]));
    }
}

void aes_decrypt_block_c(uint8_t out[16], const uint8_t in[16],
                         const aes_key *ctx)
{
    uint8_t s[16];
    int round;

    memcpy(s, in, sizeof(s));
    add_round_key(s, ctx->rk, ctx->rounds);
    for (round = ctx->rounds - 1; round > 0; round--) {
        inv_shift_rows(s);
        inv_sub_bytes(s);
        add_round_key(s, ctx->rk, round);
        inv_mix_columns(s);
    }
    inv_shift_rows(s);
    inv_sub_bytes(s);
    add_round_key(s, ctx->rk, 0);
    memcpy(out, s, sizeof(s));
}

void aes_decrypt_block(uint8_t out[16], const uint8_t in[16],
                       const aes_key *ctx)
{
    if (g_use_fast) {
        aes_decrypt_block_tt(out, in, ctx);
        return;
    }
#if defined(__riscv) && defined(AES_USE_ZVKNED)
    if (g_use_asm) {
        aes_decrypt_block_zvkned(out, in, ctx->rk_z, ctx->rounds);
        return;
    }
#endif
#if defined(__riscv) && defined(AES_USE_RVV)
    if (g_use_rvv) {
        aes_decrypt_block_rvv(out, in, ctx->rk_b, ctx->rounds);
        return;
    }
#endif
    aes_decrypt_block_c(out, in, ctx);
}

void aes_ecb_decrypt(uint8_t *out, const uint8_t *in, size_t len,
                     const aes_key *ctx)
{
    size_t off;

#if defined(__riscv) && defined(AES_USE_ZVKNED)
    if (g_use_asm && ctx->rounds == 10) {
        aes_ecb_decrypt_zvkned_128(out, in, len, ctx->rk_z);
        return;
    }
    if (g_use_asm && ctx->rounds == 12) {
        aes_ecb_decrypt_zvkned_192(out, in, len, ctx->rk_z);
        return;
    }
    if (g_use_asm && ctx->rounds == 14) {
        aes_ecb_decrypt_zvkned_256(out, in, len, ctx->rk_z);
        return;
    }
#endif
    if (g_use_fast) {
        for (off = 0; off + 16 <= len; off += 16) {
            aes_decrypt_block_tt(out + off, in + off, ctx);
        }
        return;
    }
    for (off = 0; off + 16 <= len; off += 16) {
        aes_decrypt_block(out + off, in + off, ctx);
    }
}

void aes_cbc_decrypt(uint8_t *out, const uint8_t *in, size_t len,
                     uint8_t iv[16], const aes_key *ctx)
{
    uint8_t block[16];
    uint8_t cipher[16];
    size_t off;
    int i;

#if defined(__riscv) && defined(AES_USE_ZVKNED)
    if (g_use_asm && ctx->rounds == 10) {
        aes_cbc_decrypt_zvkned_128(out, in, len, iv, ctx->rk_z);
        return;
    }
    if (g_use_asm && ctx->rounds == 12) {
        aes_cbc_decrypt_zvkned_192(out, in, len, iv, ctx->rk_z);
        return;
    }
    if (g_use_asm && ctx->rounds == 14) {
        aes_cbc_decrypt_zvkned_256(out, in, len, iv, ctx->rk_z);
        return;
    }
#endif
    if (g_use_fast) {
        for (off = 0; off + 16 <= len; off += 16) {
            uint64_t c0, c1, p0, p1, v0, v1;

            memcpy(&c0, in + off, 8);
            memcpy(&c1, in + off + 8, 8);
            memcpy(cipher, &c0, 8);
            memcpy(cipher + 8, &c1, 8);
            aes_decrypt_block_tt(block, cipher, ctx);
            memcpy(&p0, block, 8);
            memcpy(&p1, block + 8, 8);
            memcpy(&v0, iv, 8);
            memcpy(&v1, iv + 8, 8);
            p0 ^= v0;
            p1 ^= v1;
            memcpy(out + off, &p0, 8);
            memcpy(out + off + 8, &p1, 8);
            memcpy(iv, cipher, 16);
        }
        return;
    }
    for (; off + 16 <= len; off += 16) {
        memcpy(cipher, in + off, 16);
        aes_decrypt_block(block, cipher, ctx);
        for (i = 0; i < 16; i++) {
            out[off + i] = (uint8_t)(block[i] ^ iv[i]);
        }
        memcpy(iv, cipher, 16);
    }
}

void aes_cfb128_encrypt(uint8_t *out, const uint8_t *in, size_t len,
                        uint8_t iv[16], const aes_key *ctx)
{
    uint8_t block[16];
    size_t off;
    size_t i;

#if defined(__riscv) && defined(AES_USE_ZVKNED)
    {
        size_t bulk = len & ~(size_t)15;

        if (bulk != 0 && ctx->rounds == 10) {
            aes_cfb128_encrypt_zvkned_128(out, in, bulk, iv, ctx->rk_z);
        } else if (bulk != 0 && ctx->rounds == 12) {
            aes_cfb128_encrypt_zvkned_192(out, in, bulk, iv, ctx->rk_z);
        } else if (bulk != 0 && ctx->rounds == 14) {
            aes_cfb128_encrypt_zvkned_256(out, in, bulk, iv, ctx->rk_z);
        } else {
            bulk = 0;
        }
        out += bulk;
        in += bulk;
        len -= bulk;
    }
#endif
    if (g_use_fast) {
        for (off = 0; off + 16 <= len; off += 16) {
            uint64_t a, b, c, d;

            aes_encrypt_block_tt(block, iv, ctx);
            memcpy(&a, in + off, 8);
            memcpy(&b, in + off + 8, 8);
            memcpy(&c, block, 8);
            memcpy(&d, block + 8, 8);
            a ^= c;
            b ^= d;
            memcpy(out + off, &a, 8);
            memcpy(out + off + 8, &b, 8);
            memcpy(iv, out + off, 16);
        }
        if (off < len) {
            aes_encrypt_block_tt(block, iv, ctx);
            for (i = 0; off + i < len; i++) {
                out[off + i] = (uint8_t)(in[off + i] ^ block[i]);
            }
            memcpy(iv, out + off, len - off);
        }
        return;
    }
    for (off = 0; off + 16 <= len; off += 16) {
        aes_encrypt_block(block, iv, ctx);
        for (i = 0; i < 16; i++) {
            out[off + i] = (uint8_t)(in[off + i] ^ block[i]);
        }
        memcpy(iv, out + off, 16);
    }
    if (off < len) {
        aes_encrypt_block(block, iv, ctx);
        for (i = 0; off + i < len; i++) {
            out[off + i] = (uint8_t)(in[off + i] ^ block[i]);
        }
        memcpy(iv, out + off, len - off);
    }
}

void aes_cfb128_decrypt(uint8_t *out, const uint8_t *in, size_t len,
                        uint8_t iv[16], const aes_key *ctx)
{
    uint8_t block[16];
    uint8_t cipher[16];
    size_t off;
    size_t i;

#if defined(__riscv) && defined(AES_USE_ZVKNED)
    {
        size_t bulk = len & ~(size_t)15;

        if (bulk != 0 && ctx->rounds == 10) {
            aes_cfb128_decrypt_zvkned_128(out, in, bulk, iv, ctx->rk_z);
        } else if (bulk != 0 && ctx->rounds == 12) {
            aes_cfb128_decrypt_zvkned_192(out, in, bulk, iv, ctx->rk_z);
        } else if (bulk != 0 && ctx->rounds == 14) {
            aes_cfb128_decrypt_zvkned_256(out, in, bulk, iv, ctx->rk_z);
        } else {
            bulk = 0;
        }
        out += bulk;
        in += bulk;
        len -= bulk;
    }
#endif
    if (g_use_fast) {
        for (off = 0; off + 16 <= len; off += 16) {
            uint64_t a, b, c, d;

            memcpy(cipher, in + off, 16);
            aes_encrypt_block_tt(block, iv, ctx);
            memcpy(&a, cipher, 8);
            memcpy(&b, cipher + 8, 8);
            memcpy(&c, block, 8);
            memcpy(&d, block + 8, 8);
            a ^= c;
            b ^= d;
            memcpy(out + off, &a, 8);
            memcpy(out + off + 8, &b, 8);
            memcpy(iv, cipher, 16);
        }
        if (off < len) {
            aes_encrypt_block_tt(block, iv, ctx);
            for (i = 0; off + i < len; i++) {
                out[off + i] = (uint8_t)(in[off + i] ^ block[i]);
            }
        }
        return;
    }
    for (off = 0; off + 16 <= len; off += 16) {
        memcpy(cipher, in + off, 16);
        aes_encrypt_block(block, iv, ctx);
        for (i = 0; i < 16; i++) {
            out[off + i] = (uint8_t)(cipher[i] ^ block[i]);
        }
        memcpy(iv, cipher, 16);
    }
    if (off < len) {
        aes_encrypt_block(block, iv, ctx);
        for (i = 0; off + i < len; i++) {
            out[off + i] = (uint8_t)(in[off + i] ^ block[i]);
        }
    }
}

void aes_cfb8_crypt(uint8_t *out, const uint8_t *in, size_t len,
                    uint8_t iv[16], const aes_key *ctx, int encrypt)
{
    uint8_t block[16];
    uint8_t state[16];
    size_t i;
    int j;

#if defined(__riscv) && defined(AES_USE_ZVKNED)
    if (g_use_asm && ctx->rounds == 10) {
        if (encrypt) {
            aes_cfb8_crypt_zvkned_128(out, in, len, iv, ctx->rk_z);
        } else {
            aes_cfb8_crypt_zvkned_128_dec(out, in, len, iv, ctx->rk_z);
        }
        return;
    }
    if (g_use_asm && ctx->rounds == 12) {
        if (encrypt) {
            aes_cfb8_crypt_zvkned_192(out, in, len, iv, ctx->rk_z);
        } else {
            aes_cfb8_crypt_zvkned_192_dec(out, in, len, iv, ctx->rk_z);
        }
        return;
    }
    if (g_use_asm && ctx->rounds == 14) {
        if (encrypt) {
            aes_cfb8_crypt_zvkned_256(out, in, len, iv, ctx->rk_z);
        } else {
            aes_cfb8_crypt_zvkned_256_dec(out, in, len, iv, ctx->rk_z);
        }
        return;
    }
#endif
    memcpy(state, iv, 16);
    for (i = 0; i < len; i++) {
        aes_encrypt_block(block, state, ctx);
        out[i] = (uint8_t)(in[i] ^ block[0]);
        for (j = 0; j < 15; j++) {
            state[j] = state[j + 1];
        }
        state[15] = encrypt ? out[i] : in[i];
    }
    memcpy(iv, state, 16);
}

void aes_cfb1_crypt(uint8_t *out_bits, const uint8_t *in_bits, size_t nbits,
                    uint8_t iv[16], const aes_key *ctx, int encrypt)
{
    uint8_t block[16];
    uint8_t state[16];
    size_t i;
    int j;

    memcpy(state, iv, 16);
    for (i = 0; i < nbits; i++) {
        uint8_t feedback;
        aes_encrypt_block(block, state, ctx);
        out_bits[i] = (uint8_t)(in_bits[i] ^ ((block[0] >> 7) & 1));
        feedback = encrypt ? out_bits[i] : in_bits[i];
        for (j = 0; j < 15; j++) {
            state[j] = (uint8_t)((state[j] << 1) | (state[j + 1] >> 7));
        }
        state[15] = (uint8_t)((state[15] << 1) | feedback);
    }
    memcpy(iv, state, 16);
}

void aes_ofb_crypt(uint8_t *out, const uint8_t *in, size_t len,
                   uint8_t iv[16], const aes_key *ctx)
{
    uint8_t block[16];
    size_t off;
    size_t i;

#if defined(__riscv) && defined(AES_USE_ZVKNED)
    {
        size_t bulk = len & ~(size_t)15;

        if (bulk != 0 && ctx->rounds == 10) {
            aes_ofb_crypt_zvkned_128(out, in, bulk, iv, ctx->rk_z);
        } else if (bulk != 0 && ctx->rounds == 12) {
            aes_ofb_crypt_zvkned_192(out, in, bulk, iv, ctx->rk_z);
        } else if (bulk != 0 && ctx->rounds == 14) {
            aes_ofb_crypt_zvkned_256(out, in, bulk, iv, ctx->rk_z);
        } else {
            bulk = 0;
        }
        out += bulk;
        in += bulk;
        len -= bulk;
    }
#endif
    for (off = 0; off < len; off += 16) {
        size_t n = len - off < 16 ? len - off : 16;
        aes_encrypt_block(block, iv, ctx);
        memcpy(iv, block, 16);
        for (i = 0; i < n; i++) {
            out[off + i] = (uint8_t)(in[off + i] ^ block[i]);
        }
    }
}
